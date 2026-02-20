# 🧬 AKT1 Inhibitor RNA-seq (MDS) — Snakemake + Salmon + DESeq2

A reproducible RNA-seq pipeline to evaluate how an **AKT1 inhibitor (PU-001)** affects **immune checkpoint gene expression** in **myelodysplastic syndrome (MDS)**-related cell line models.

Built for: **fast QC → trimming → Salmon quantification → tximport → DESeq2 → checkpoint summary → MultiQC**.

---

## ✨ Project goals

✅ Quantify transcript/gene expression across conditions  
✅ Test differential expression (CTRL vs PU-001)  
✅ Produce interpretable outputs for a thesis:
- PCA / sample QC
- DE tables per cell line
- A focused **checkpoint gene panel summary** (e.g., *CD274/PD-L1, HAVCR2/TIM-3*, etc.)

---

## 🧪 Experimental design (example)

- **Cell lines:** KG-1, Mono-Mac-1, THP-1  
- **Conditions:** CTRL vs PU-001  
- **Replicates:** 3 per condition  
- **Total samples:** 18 paired-end FASTQ files

---

## 🚀 Quick start

### Prepare input FASTQ files
Place paired-end FASTQ in:
```text
data/raw/
  SAMPLE_1.fq.gz
  SAMPLE_2.fq.gz
```

Sample names must match the `sample` column in `config/samples.tsv`.

---

## ⚙️ Configuration

### `config/samples.tsv`
A minimal sample sheet:

```tsv
sample	cell_line	treatment	r1	r2
T17	THP-1	CTRL	data/raw/T17_1.fq.gz	data/raw/T17_2.fq.gz
T18	THP-1	PU001	data/raw/T18_1.fq.gz	data/raw/T18_2.fq.gz
...
```

### `config/config.yaml`
Controls key paths + parameters (threads, reference files, etc.).

---

## 🧰 Install dependencies (recommended: Snakemake with conda)

### Option A — Let Snakemake manage environments (best reproducibility)
Install Snakemake once:
```bash
mamba create -n snakemake -c conda-forge -c bioconda snakemake
conda activate snakemake
```

Then run with:
```bash
snakemake -s workflow/Snakefile --use-conda --cores 4 -p --latency-wait 60
```

### Option B — Use your own environment
If you already have an environment (e.g. `rnaseq`), ensure it contains:
- `fastp`, `fastqc`, `salmon`, `multiqc`
- R + packages: `tximport`, `DESeq2`, `pheatmap`, `ggplot2`, etc.

---

## 🐍 Pipeline overview (Salmon-based)

### ✅ Step-by-step execution order

1. **QC (raw reads)** 🧾  
   Rule: `fastqc_raw`  
   Output: `results/fastqc_raw/`

2. **Trimming** ✂️  
   Rule: `fastp_trim`  
   Output: `results/trimmed/{sample}_1.fq.gz`, `results/trimmed/{sample}_2.fq.gz`

3. **Salmon index** 🧱  
   Rule: `salmon_index`  
   Output: `ref/salmon_index_gencode_v49/` (or configured path)

4. **Salmon quant (per sample)** ⚡  
   Rule: `salmon_quant`  
   Output:
   - `results/salmon/{sample}/quant.sf`
   - `results/salmon/{sample}/salmon.log`

5. **tximport + DESeq2 (gene-level analysis)** 📊  
   Script: `workflow/scripts/run_deseq2_tximport.R`  
   Outputs (examples):
   - `results/deseq2/deseq2_results_all.tsv`
   - `results/deseq2/by_cell_line/<cell_line>_CTRL_vs_PU001.tsv`
   - `results/deseq2/pca_plot.png`
   - `results/deseq2/checkpoints_summary.tsv`

6. **MultiQC report** 📦  
   Rule: `multiqc`  
   Output: `results/multiqc/multiqc_report.html`

---

## 🧾 How to run (common commands)

### Run everything
```bash
snakemake -s workflow/Snakefile --cores 4 -p --latency-wait 60
```

### Force rerun incomplete jobs (after interruption)
```bash
snakemake -s workflow/Snakefile --cores 4 -p --rerun-incomplete --latency-wait 60
```

### Mark manually-provided outputs as complete (cleanup metadata)
If you generated `results/salmon/<sample>/quant.sf` externally (e.g., Colab) and Snakemake says “incomplete”:
```bash
snakemake --cleanup-metadata results/salmon/T17/quant.sf results/salmon/T17/salmon.log
```

---

## 🧠 Checkpoint gene panel

Checkpoint genes are listed in:
- `config/checkpoint_genes.txt`

The pipeline generates:
- `results/deseq2/checkpoints_summary.tsv`

This table is designed to quickly answer:
➡️ “Does PU-001 reduce checkpoint gene expression?”

---

## 👤 Author

Maintained by: **Alexey Demaev**  
Affiliation: **ITMO University**  
Contact: **<alexeydemq@gmail.com>**
