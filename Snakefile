# -*- coding: utf-8 -*-
# Salmon-based RNA-seq pipeline: FastQC -> fastp -> Salmon -> tximport+DESeq2 -> MultiQC
# Assumes paired-end reads.
#
# Required inputs:
#   - metadata/samples.tsv  (columns: sample_id, fq1, fq2, cell_line, treatment, replicate)
#   - ref/ transcriptome fasta + GTF (GENCODE release should match)
#
# Outputs:
#   - results/salmon/{sample}/quant.sf
#   - results/deseq2/* (DE tables, PCA, checkpoint summary, heatmap)
#   - results/multiqc/multiqc_report.html

import csv
from pathlib import Path

configfile: "config.yaml"

SAMPLES_TSV = config["samples_tsv"]

def load_samples(tsv_path):
    rows = []
    with open(tsv_path, "r", newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for r in reader:
            rows.append(r)
    need = {"sample_id","fq1","fq2","cell_line","treatment","replicate"}
    if not rows:
        raise ValueError(f"samples.tsv seems empty: {tsv_path}")
    missing = need - set(rows[0].keys())
    if missing:
        raise ValueError(f"samples.tsv missing columns: {sorted(missing)}")
    return rows

_rows = load_samples(SAMPLES_TSV)
SAMPLES = [r["sample_id"] for r in _rows]
CELL_LINES = sorted(set(r["cell_line"] for r in _rows))

fq1_map = {r["sample_id"]: r["fq1"] for r in _rows}
fq2_map = {r["sample_id"]: r["fq2"] for r in _rows}

def fq1(wc): return fq1_map[wc.sample]
def fq2(wc): return fq2_map[wc.sample]

# ---------- Paths from config ----------
REF_GTF = config["ref"]["gtf"]
TRANSCRIPTS_FA = config["ref"]["transcriptome_fasta"]
SALMON_INDEX = config["ref"]["salmon_index"]

THREADS_FASTQC = int(config["threads"].get("fastqc", 4))
THREADS_FASTP  = int(config["threads"].get("fastp", 4))
THREADS_SALMON = int(config["threads"].get("salmon", 4))

FASTP_EXTRA = config.get("fastp", {}).get("extra", "--detect_adapter_for_pe")
SALMON_EXTRA = config.get("salmon", {}).get("extra", "--gcBias --seqBias --posBias")
SALMON_KMER  = int(config.get("salmon", {}).get("kmer", 31))

CHECKPOINT_GENES = config.get("checkpoint_genes", [
    "CD274","PDCD1LG2","PDCD1","CTLA4","LAG3","TIGIT","HAVCR2","BTLA","VSIR","IDO1","CD276","VTCN1","CD80","CD86"
])

# ---------- Final targets ----------
rule all:
    input:
        "results/multiqc/multiqc_report.html",
        "results/deseq2/checkpoints_summary.tsv",
        "results/deseq2/DE_all_cell_lines.tsv",
        expand("results/deseq2/DE_{cl}.tsv", cl=CELL_LINES),
        expand("results/salmon/{sample}/quant.sf", sample=SAMPLES)

# ---------- QC ----------
rule qc_raw:
    input:
        expand("results/fastqc/raw/{sample}_1_fastqc.html", sample=SAMPLES),
        expand("results/fastqc/raw/{sample}_2_fastqc.html", sample=SAMPLES)

rule fastqc_raw:
    input:
        r1=fq1,
        r2=fq2
    output:
        html1="results/fastqc/raw/{sample}_1_fastqc.html",
        zip1 ="results/fastqc/raw/{sample}_1_fastqc.zip",
        html2="results/fastqc/raw/{sample}_2_fastqc.html",
        zip2 ="results/fastqc/raw/{sample}_2_fastqc.zip",
    threads: THREADS_FASTQC
    shell:
        r"""
        mkdir -p results/fastqc/raw
        fastqc -t {threads} -o results/fastqc/raw {input.r1} {input.r2}

        b1=$(basename "{input.r1}")
        b2=$(basename "{input.r2}")
        b1=${{b1%.fq.gz}}; b1=${{b1%.fastq.gz}}
        b2=${{b2%.fq.gz}}; b2=${{b2%.fastq.gz}}

        mv -f "results/fastqc/raw/${{b1}}_fastqc.html" "{output.html1}"
        mv -f "results/fastqc/raw/${{b1}}_fastqc.zip"  "{output.zip1}"
        mv -f "results/fastqc/raw/${{b2}}_fastqc.html" "{output.html2}"
        mv -f "results/fastqc/raw/${{b2}}_fastqc.zip"  "{output.zip2}"
        """

# ---------- Trimming ----------
rule fastp_trim:
    input:
        r1=fq1,
        r2=fq2
    output:
        r1="results/trimmed/{sample}_1.fq.gz",
        r2="results/trimmed/{sample}_2.fq.gz",
        html="results/fastp/{sample}.html",
        json="results/fastp/{sample}.json"
    threads: THREADS_FASTP
    shell:
        r"""
        mkdir -p results/trimmed results/fastp
        fastp -w {threads} \
          -i {input.r1} -I {input.r2} \
          -o {output.r1} -O {output.r2} \
          -h {output.html} -j {output.json} \
          {FASTP_EXTRA}
        """

# ---------- Salmon index ----------
rule salmon_index:
    input:
        fa=TRANSCRIPTS_FA
    output:
        directory(SALMON_INDEX)
    threads: 4
    shell:
        r"""
        mkdir -p {output}
        salmon index -t {input.fa} -i {output} -k {SALMON_KMER}
        """

# ---------- Salmon quant ----------
rule salmon_quant:
    input:
        idx=SALMON_INDEX,
        r1="results/trimmed/{sample}_1.fq.gz",
        r2="results/trimmed/{sample}_2.fq.gz"
    output:
        quant="results/salmon/{sample}/quant.sf",
        log="results/salmon/{sample}/salmon.log"
    threads: THREADS_SALMON
    shell:
        r"""
        mkdir -p results/salmon/{wildcards.sample}
        salmon quant -i {input.idx} -l A \
          -1 {input.r1} -2 {input.r2} \
          -p {threads} \
          {SALMON_EXTRA} \
          -o results/salmon/{wildcards.sample} \
          &> {output.log}
        test -s {output.quant}
        """

# ---------- tximport + DESeq2 + checkpoint hypothesis check ----------
rule deseq2_tximport:
    input:
        samples=SAMPLES_TSV,
        gtf=REF_GTF,
        quants=expand("results/salmon/{sample}/quant.sf", sample=SAMPLES)
    output:
        "results/deseq2/checkpoints_summary.tsv",
        "results/deseq2/DE_all_cell_lines.tsv",
        expand("results/deseq2/DE_{cl}.tsv", cl=CELL_LINES),
        "results/deseq2/pca_vst_all.png",
        "results/deseq2/checkpoints_heatmap.png"
    params:
        quant_dir="results/salmon",
        outdir="results/deseq2",
        checkpoints=",".join(CHECKPOINT_GENES)
    shell:
        r"""
        mkdir -p {params.outdir}
        Rscript scripts/tximport_deseq2.R \
          --samples {input.samples} \
          --gtf {input.gtf} \
          --quant_dir {params.quant_dir} \
          --outdir {params.outdir} \
          --checkpoints "{params.checkpoints}"
        """

# ---------- MultiQC ----------
rule multiqc:
    input:
        "results/deseq2/checkpoints_summary.tsv"
    output:
        "results/multiqc/multiqc_report.html"
    threads: 1
    shell:
        r"""
        mkdir -p results/multiqc
        multiqc -o results/multiqc results
        """
