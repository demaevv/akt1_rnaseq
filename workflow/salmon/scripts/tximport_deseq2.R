#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(tximport)
  library(DESeq2)
  library(ggplot2)
  library(pheatmap)
})

opt_list <- list(
  make_option(c("--samples"), type="character", help="metadata/samples.tsv"),
  make_option(c("--gtf"), type="character", help="GENCODE .gtf (same release as transcriptome)"),
  make_option(c("--quant_dir"), type="character", default="results/salmon", help="Base dir with per-sample Salmon folders"),
  make_option(c("--outdir"), type="character", default="results/deseq2", help="Output directory"),
  make_option(c("--checkpoints"), type="character", default="", help="Comma-separated gene symbols to test")
)
opt <- parse_args(OptionParser(option_list=opt_list))

if (is.null(opt$samples) || is.null(opt$gtf)) stop("Required: --samples and --gtf")
dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)

# -----------------------------
# Load sample sheet
# -----------------------------
samp <- fread(opt$samples, sep="\t", header=TRUE)
need_cols <- c("sample_id","cell_line","treatment","replicate")
miss <- setdiff(need_cols, names(samp))
if (length(miss) > 0) stop(paste("samples.tsv missing columns:", paste(miss, collapse=",")))

# Build quant.sf paths
samp[, quant_path := file.path(opt$quant_dir, sample_id, "quant.sf")]
if (any(!file.exists(samp$quant_path))) {
  missing <- samp[!file.exists(quant_path), .(sample_id, quant_path)]
  fwrite(missing, file.path(opt$outdir, "missing_quant_paths.tsv"), sep="\t")
  stop("Some quant.sf files are missing. See results/deseq2/missing_quant_paths.tsv")
}
files <- samp$quant_path
names(files) <- samp$sample_id

# -----------------------------
# Build tx2gene + gene_name mapping from GTF
# -----------------------------
gtf_cmd <- sprintf("grep -v '^#' %s", shQuote(opt$gtf))
gtf <- data.table::fread(cmd = gtf_cmd, sep = "\t", header = FALSE, quote = "", fill = TRUE)

if (ncol(gtf) < 9) stop("GTF parse failed: expected >= 9 columns")
if (ncol(gtf) > 9) gtf <- gtf[, 1:9]

setnames(gtf, c("seqname","source","feature","start","end","score","strand","frame","attribute"))

extract_attr <- function(x, key) {
  pat <- paste0(key, ' "([^"]+)"')
  ifelse(grepl(pat, x), sub(paste0(".*", pat, ".*"), "\\1", x), NA_character_)
}

gtf[, transcript_id := extract_attr(attribute, "transcript_id")]
gtf[, gene_id       := extract_attr(attribute, "gene_id")]
gtf[, gene_name     := extract_attr(attribute, "gene_name")]

tx2gene <- unique(gtf[!is.na(transcript_id) & !is.na(gene_id), .(TXNAME=transcript_id, GENEID=gene_id)])
gene_annot <- unique(gtf[!is.na(gene_id) & !is.na(gene_name), .(GENEID=gene_id, gene_name)])

strip_ver <- function(x) sub("\\..*$", "", x)

tx2gene[, TXNAME := strip_ver(TXNAME)]
tx2gene[, GENEID := strip_ver(GENEID)]
tx2gene <- unique(tx2gene)

gene_annot[, GENEID := strip_ver(GENEID)]
setorder(gene_annot, GENEID)
gene_annot <- gene_annot[, .(gene_name = gene_name[1]), by = GENEID]

fwrite(tx2gene, file.path(opt$outdir, "tx2gene.tsv"), sep="\t")
fwrite(gene_annot, file.path(opt$outdir, "gene_annot.tsv"), sep="\t")

# -----------------------------
# tximport -> gene-level counts
# -----------------------------
txi <- tximport(
  files,
  type="salmon",
  tx2gene=tx2gene,
  ignoreTxVersion=TRUE,
  countsFromAbundance="lengthScaledTPM",
  dropInfReps=TRUE
)

# -----------------------------
# colData
# -----------------------------
coldata <- as.data.frame(samp[, .(sample_id, cell_line, treatment, replicate)])
rownames(coldata) <- coldata$sample_id
coldata$cell_line <- factor(coldata$cell_line)
coldata$treatment <- factor(coldata$treatment)
if ("CTRL" %in% levels(coldata$treatment)) coldata$treatment <- relevel(coldata$treatment, "CTRL")

coldata <- coldata[colnames(txi$counts), , drop=FALSE]
stopifnot(identical(colnames(txi$counts), rownames(coldata)))

treated <- setdiff(levels(coldata$treatment), "CTRL")
if (length(treated) != 1) stop("Expected exactly one treated level besides CTRL. Found: ", paste(treated, collapse=","))
treated <- treated[1]

write_de <- function(res, out_path) {
  dt <- as.data.table(res, keep.rownames = "GENEID")
  dt <- merge(dt, gene_annot, by="GENEID", all.x=TRUE)
  setcolorder(dt, c("GENEID","gene_name", setdiff(names(dt), c("GENEID","gene_name"))))
  fwrite(dt[order(padj)], out_path, sep="\t")
}

# Helper: subset txi to the same samples as colData
subset_txi <- function(txi, sample_ids) {
  txi2 <- txi
  txi2$abundance <- txi$abundance[, sample_ids, drop=FALSE]
  txi2$counts    <- txi$counts[,    sample_ids, drop=FALSE]
  txi2$length    <- txi$length[,    sample_ids, drop=FALSE]
  txi2
}

# -----------------------------
# DESeq2 (all cell lines; controls for cell_line)
# -----------------------------
dds_all <- DESeqDataSetFromTximport(txi, colData=coldata, design=~ cell_line + treatment)
dds_all <- dds_all[rowSums(counts(dds_all)) >= 10, ]
dds_all <- DESeq(dds_all)
res_all <- results(dds_all, contrast=c("treatment", treated, "CTRL"))
write_de(res_all, file.path(opt$outdir, "DE_all_cell_lines.tsv"))

# PCA (VST)
vsd <- vst(dds_all, blind=TRUE)
pca <- plotPCA(vsd, intgroup=c("cell_line","treatment"), returnData=TRUE)
percentVar <- round(100 * attr(pca, "percentVar"))
p <- ggplot(pca, aes(PC1, PC2, color=cell_line, shape=treatment)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ", percentVar[1], "%")) +
  ylab(paste0("PC2: ", percentVar[2], "%")) +
  theme_bw()
ggsave(filename=file.path(opt$outdir, "pca_vst_all.png"), plot=p, width=7, height=5, dpi=200)

# -----------------------------
# Per-cell-line DE (treatment effect)
# -----------------------------
for (cl in levels(coldata$cell_line)) {
  idx <- coldata$cell_line == cl
  sample_ids <- rownames(coldata)[idx]

  txi_cl <- subset_txi(txi, sample_ids)
  cd_cl  <- coldata[sample_ids, , drop=FALSE]

  dds <- DESeqDataSetFromTximport(txi_cl, colData=cd_cl, design=~ treatment)
  dds <- dds[rowSums(counts(dds)) >= 10, ]
  dds <- DESeq(dds)
  res <- results(dds, contrast=c("treatment", treated, "CTRL"))
  write_de(res, file.path(opt$outdir, paste0("DE_", cl, ".tsv")))
}

# -----------------------------
# Hypothesis check: checkpoint genes downregulated after AKT1 inhibition
# -----------------------------
ck_syms <- character(0)
if (nchar(opt$checkpoints) > 0) {
  ck_syms <- trimws(unlist(strsplit(opt$checkpoints, ",")))
  ck_syms <- ck_syms[nzchar(ck_syms)]
}

ck_geneids <- gene_annot[gene_name %in% ck_syms, unique(GENEID)]
missing_syms <- setdiff(ck_syms, gene_annot$gene_name)
if (length(missing_syms) > 0) {
  fwrite(data.table(missing_symbol=missing_syms),
         file.path(opt$outdir, "checkpoints_missing_symbols.tsv"), sep="\t")
}

summ_list <- list()
for (cl in levels(coldata$cell_line)) {
  idx <- coldata$cell_line == cl
  sample_ids <- rownames(coldata)[idx]

  txi_cl <- subset_txi(txi, sample_ids)
  cd_cl  <- coldata[sample_ids, , drop=FALSE]

  dds <- DESeqDataSetFromTximport(txi_cl, colData=cd_cl, design=~ treatment)
  dds <- dds[rowSums(counts(dds)) >= 10, ]
  dds <- DESeq(dds)
  res <- results(dds, contrast=c("treatment", treated, "CTRL"))

  dt <- as.data.table(res, keep.rownames="GENEID")
  dt <- merge(dt, gene_annot, by="GENEID", all.x=TRUE)
  dt[, cell_line := cl]
  dt <- dt[GENEID %in% ck_geneids]
  summ_list[[cl]] <- dt
}
ck_sum <- rbindlist(summ_list, fill=TRUE)
setcolorder(ck_sum, c("cell_line","GENEID","gene_name","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj"))
fwrite(ck_sum[order(cell_line, padj)], file.path(opt$outdir, "checkpoints_summary.tsv"), sep="\t")

# Heatmap (VST)
if (length(ck_geneids) > 1) {
  mat <- assay(vsd)

  mat_ids_clean <- sub("\\..*$", "", rownames(mat))
  ck_ids_clean  <- sub("\\..*$", "", ck_geneids)

  keep_clean <- intersect(ck_ids_clean, mat_ids_clean)

  if (length(keep_clean) > 1) {
    idx <- mat_ids_clean %in% keep_clean
    hm  <- mat[idx, , drop = FALSE]

    rownames(hm) <- mat_ids_clean[idx]

    # Column annotations
    ann_col <- data.frame(
      cell_line = coldata$cell_line,
      treatment = coldata$treatment
    )
    rownames(ann_col) <- rownames(coldata)

    # ---- gene_name labels for rows ----
    labels_row <- rownames(hm)  # fallback: gene_id

    if (exists("gene_annot")) {
      ga <- as.data.frame(gene_annot)

      id_candidates   <- c("gene_id","GENEID","GeneID","geneid","geneId","gid")
      name_candidates <- c("gene_name","GENENAME","GeneName","genename","symbol","SYMBOL","hgnc_symbol","gene")

      id_col   <- id_candidates[id_candidates %in% names(ga)][1]
      name_col <- name_candidates[name_candidates %in% names(ga)][1]

      if (!is.na(id_col) && !is.na(name_col)) {
        ids   <- sub("\\..*$", "", as.character(ga[[id_col]]))
        names <- as.character(ga[[name_col]])

        # take first non-empty gene_name per gene_id
        tmp <- stats::aggregate(
          names,
          by = list(gene_id = ids),
          FUN = function(x) x[which(nzchar(x))[1]]
        )
        colnames(tmp)[2] <- "gene_name"

        mapped <- tmp$gene_name[match(rownames(hm), tmp$gene_id)]
        mapped[is.na(mapped) | mapped == ""] <- rownames(hm)[is.na(mapped) | mapped == ""]
        labels_row <- make.unique(mapped)
      }
    }
    # ------------------------------------------------------------

    png(file.path(opt$outdir, "checkpoints_heatmap.png"), width = 1200, height = 900)
    pheatmap(
      hm,
      scale = "row",
      annotation_col = ann_col,
      show_colnames = FALSE,
      fontsize_row = 9,
      labels_row = labels_row,
      main = paste0("Checkpoint genes (VST; scaled) — ", treated, " vs CTRL")
    )
    dev.off()
  }
}

cat("OK\n")
