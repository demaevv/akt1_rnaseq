# 🧬 AKT1 RNA-seq Project  
## Molecular mechanisms of PU-001 response in MDS-relevant myeloid cell models

This repository contains a reproducible RNA-seq analysis framework for studying the effects of the novel AKT1 inhibitor **PU-001** on immune checkpoint-associated transcriptional and pathway-level programs in myelodysplastic syndrome (MDS)-relevant myeloid cell models.

The project integrates two independent RNA-seq quantification strategies:

- ⭐ **STAR + featureCounts + DESeq2**
- 🐟 **Salmon + tximport + DESeq2**

and extends them with downstream:

- 🧠 checkpoint-focused interpretation
- 🧭 pathway enrichment analysis
- 📊 ORA / Enrichr / preranked GSEA
- 🔬 mechanistic hypothesis generation

---

## 🎯 Biological question

Does pharmacological AKT1 inhibition by **PU-001** remodel immune checkpoint-associated expression programs in MDS-relevant myeloid cell models?

More specifically, the project asks whether PU-001 affects:

- **PD-L1 axis** — `CD274`, `PDCD1LG2`
- **TIM-3 axis** — `HAVCR2`, `LGALS9`
- **B7 family molecules** — `CD80`, `CD86`, `CD276`
- **alternative checkpoint regulators** — `VSIR`, `CTLA4`, `LAG3`, `IDO1`, `BTLA`
- pathway-level programs related to **PI3K/AKT/mTOR**, **IFN/JAK-STAT**, **NF-κB**, antigen presentation, cell cycle, translation, and stress adaptation

---

## 🧪 Experimental design

The RNA-seq dataset contains **18 paired-end libraries**:

| Cell line | Control | PU-001, 200 μM, 24 h | Replicates |
|---|---|---|---|
| KG-1 | K1, K3, K5 | K2, K4, K6 | 3 + 3 |
| Mono-Mac-1 | M7, M9, M11 | M8, M10, M12 | 3 + 3 |
| THP-1 | T13, T15, T17 | T14, T16, T18 | 3 + 3 |

The cell lines represent distinct myeloid / monocytic differentiation contexts:

- **KG-1** — immature AML-like model with elevated basal TIM-3 / PD-L1
- **Mono-Mac-1** — myelomonocytic model
- **THP-1** — monocytic model

For THP-1 and Mono-Mac-1, IFN-γ was used to induce checkpoint expression before evaluating the effect of PU-001.

---

## ⚙️ Workflow overview

### 1️⃣ Input data

Place paired-end FASTQ files into:

```text
data/raw/
```

Example:

```text
data/raw/Unknown_CQ888-001U0001_1.fq.gz
data/raw/Unknown_CQ888-001U0001_2.fq.gz
```

Sample metadata are defined in:

```text
config/samples.tsv
```

Recommended shared sample-sheet format:

```tsv
sample_id	client_id	cell_line	treatment	replicate	fq1	fq2
K1	CQ888-001U0001	KG-1	CTRL	1	data/raw/Unknown_CQ888-001U0001_1.fq.gz	data/raw/Unknown_CQ888-001U0001_2.fq.gz
K2	CQ888-001U0002	KG-1	PU001	1	data/raw/Unknown_CQ888-001U0002_1.fq.gz	data/raw/Unknown_CQ888-001U0002_2.fq.gz
K3	CQ888-001U0003	KG-1	CTRL	2	data/raw/Unknown_CQ888-001U0003_1.fq.gz	data/raw/Unknown_CQ888-001U0003_2.fq.gz
K4	CQ888-001U0004	KG-1	PU001	2	data/raw/Unknown_CQ888-001U0004_1.fq.gz	data/raw/Unknown_CQ888-001U0004_2.fq.gz
```

---

## ⭐ STAR-based branch

The STAR branch implements a classical genome-alignment workflow:

```text
FASTQ
  ↓
FastQC
  ↓
fastp
  ↓
STAR genome alignment
  ↓
samtools sort/index
  ↓
featureCounts
  ↓
DESeq2
  ↓
checkpoint summary + heatmaps + PCA
```

### Main files

```text
workflow/star/Snakefile
workflow/star/scripts/get_ref.sh
workflow/star/scripts/featurecounts_to_tsv.py
workflow/star/scripts/deseq2_pipeline.R
config/star.yaml
```

### Run STAR workflow

```bash
snakemake \
  -s workflow/star/Snakefile \
  --configfile config/star.yaml \
  --cores 12 \
  -p \
  --rerun-incomplete \
  --latency-wait 60
```

### Expected outputs

```text
results/star/fastqc/
results/star/trimmed/
results/star/bam/
results/star/counts/gene_counts.tsv
results/star/deseq2/DE_KG-1_PU001_vs_CTRL.tsv
results/star/deseq2/DE_Mono-Mac-1_PU001_vs_CTRL.tsv
results/star/deseq2/DE_THP-1_PU001_vs_CTRL.tsv
results/star/deseq2/checkpoints_summary.tsv
results/star/multiqc/multiqc_report.html
```

---

## 🐟 Salmon-based branch

The Salmon branch implements transcript-aware quantification:

```text
FASTQ
  ↓
FastQC
  ↓
fastp
  ↓
Salmon transcript quantification
  ↓
tximport
  ↓
DESeq2
  ↓
checkpoint summary + PCA
```

### Main files

```text
workflow/salmon/Snakefile
workflow/salmon/scripts/get_transcriptome.sh
workflow/salmon/scripts/tximport_deseq2.R
config/salmon.yaml
```

### Run Salmon workflow

```bash
snakemake \
  -s workflow/salmon/Snakefile \
  --configfile config/salmon.yaml \
  --cores 12 \
  -p \
  --rerun-incomplete \
  --latency-wait 60
```

### Expected outputs

```text
results/salmon/fastqc/
results/salmon/trimmed/
results/salmon/quant/
results/salmon/deseq2/DE_KG-1_PU001_vs_CTRL.tsv
results/salmon/deseq2/DE_Mono-Mac-1_PU001_vs_CTRL.tsv
results/salmon/deseq2/DE_THP-1_PU001_vs_CTRL.tsv
results/salmon/deseq2/checkpoints_summary.tsv
results/salmon/multiqc/multiqc_report.html
```

---

## 🧠 Checkpoint gene panel

The checkpoint gene panel is stored in:

```text
config/checkpoint_genes.txt
```

Core genes:

```text
CD274
HAVCR2
LAG3
CTLA4
CD276
VSIR
CD80
CD86
PDCD1
PDCD1LG2
IDO1
BTLA
LGALS9
TNFRSF14
```

The checkpoint summary is designed to answer:

> Does PU-001 reduce checkpoint-associated gene expression, and is the effect reproducible across quantification strategies?

---

## 📊 Cross-pipeline comparison

The key analytical principle of this project is **cross-pipeline reproducibility**.

A checkpoint change is treated as high-confidence if:

1. the gene is detected in both STAR and Salmon branches;
2. the direction of change is the same;
3. the adjusted p-value supports statistical significance;
4. the result is biologically interpretable in the context of pathway analysis.

---

## 🧭 Pathway enrichment analysis

Pathway analysis is performed after obtaining DESeq2 tables.

Scripts:

```text
workflow/pathway_enrichment/scripts/pathway_analysis_kg1.py
workflow/pathway_enrichment/scripts/pathway_analysis_monomac1.py
workflow/pathway_enrichment/scripts/pathway_analysis_thp1.py
```

Each script performs:

- filtering by `padj < 0.05` and `|log2FC| ≥ 1`
- top-200 UP and DOWN gene selection
- ORA using MSigDB-style collections
- Enrichr-based cross-validation
- preranked GSEA using the full ranked gene list
- dotplots, barplots, GSEA running plots
- markdown report generation

### Recommended input placement

The pathway scripts should read STAR-based DESeq2 tables from:

```text
results/star/deseq2/
```

Recommended paths:

```text
results/star/deseq2/DE_KG-1_PU001_vs_CTRL.tsv
results/star/deseq2/DE_Mono-Mac-1_PU001_vs_CTRL.tsv
results/star/deseq2/DE_THP-1_PU001_vs_CTRL.tsv
```

### Run pathway analysis

```bash
python workflow/pathway_enrichment/scripts/pathway_analysis_kg1.py
python workflow/pathway_enrichment/scripts/pathway_analysis_monomac1.py
python workflow/pathway_enrichment/scripts/pathway_analysis_thp1.py
```

Expected output directories:

```text
results/pathway_enrichment/KG1/
results/pathway_enrichment/MonoMac1/
results/pathway_enrichment/THP1/
```

---

## 🚀 Quick start

### 1. Clone the repository

```bash
git clone https://github.com/demaevv/akt1_rnaseq.git
cd akt1_rnaseq
```

### 2. Create the main environment

```bash
mamba env create -f environment.yml
conda activate akt1_rnaseq
```

### 3. Prepare input data

Place FASTQ files into:

```text
data/raw/
```

Update:

```text
config/samples.tsv
```

### 4. Run STAR branch

```bash
snakemake -s workflow/star/Snakefile --configfile config/star.yaml --cores 12 -p
```

### 5. Run Salmon branch

```bash
snakemake -s workflow/salmon/Snakefile --configfile config/salmon.yaml --cores 12 -p
```

### 6. Run pathway enrichment

```bash
python workflow/pathway_enrichment/scripts/pathway_analysis_kg1.py
python workflow/pathway_enrichment/scripts/pathway_analysis_monomac1.py
python workflow/pathway_enrichment/scripts/pathway_analysis_thp1.py
```

---

## 👤 Author

**Alexey Demaev**  
M4235, Bioinformatics and Systems Biology  
ITMO University

GitHub: [@demaevv](https://github.com/demaevv)
