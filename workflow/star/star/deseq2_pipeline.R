#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
counts_tsv <- args[1]
meta_tsv   <- args[2]
gtf_path   <- args[3]
outdir     <- args[4]

suppressPackageStartupMessages({
  library(DESeq2)
  library(data.table)
  library(ggplot2)
  library(pheatmap)
  library(rtracklayer)
})

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(outdir, "figures"), showWarnings = FALSE, recursive = TRUE)

# -----------------------------
# Load data
# -----------------------------
counts <- fread(counts_tsv)  # gene_id + sample columns
meta   <- fread(meta_tsv)

stopifnot("sample_id" %in% colnames(meta))
stopifnot(all(meta$sample_id %in% colnames(counts)))

meta <- as.data.frame(meta)
meta$cell_line  <- factor(meta$cell_line)
meta$treatment  <- factor(meta$treatment, levels = c("CTRL","PU001"))
meta$replicate  <- factor(meta$replicate)

rownames(meta) <- meta$sample_id

# Counts matrix
count_mat <- as.matrix(counts[, -1])
rownames(count_mat) <- counts$gene_id
count_mat <- count_mat[, meta$sample_id, drop=FALSE]

# -----------------------------
# Map gene_id -> gene_name from GTF
# -----------------------------
gtf <- import(gtf_path)
gtf_genes <- unique(mcols(gtf)[, c("gene_id","gene_name")])
gtf_genes <- gtf_genes[!is.na(gtf_genes$gene_id) & !is.na(gtf_genes$gene_name),]
gene_map <- setNames(as.character(gtf_genes$gene_name), as.character(gtf_genes$gene_id))

gene_name <- gene_map[rownames(count_mat)]
gene_name[is.na(gene_name)] <- rownames(count_mat)
names(gene_name) <- rownames(count_mat)

# -----------------------------
# DESeq2: global model (controls cell_line)
# -----------------------------
dds_all <- DESeqDataSetFromMatrix(
  countData = count_mat,
  colData   = meta,
  design    = ~ cell_line + treatment
)

keep <- rowSums(counts(dds_all) >= 10) >= 3
dds_all <- dds_all[keep,]
dds_all <- DESeq(dds_all)

res_all <- results(dds_all, contrast = c("treatment","PU001","CTRL"))
res_all <- as.data.frame(res_all)
res_all$gene_id <- rownames(res_all)
res_all$gene_name <- gene_name[res_all$gene_id]
res_all <- res_all[order(res_all$padj),]
write.table(res_all, file=file.path(outdir,"DE_allCellLines_treatment.tsv"),
            sep="\t", quote=FALSE, row.names=FALSE)

# PCA
vsd <- vst(dds_all, blind=FALSE)
pca <- plotPCA(vsd, intgroup=c("cell_line","treatment"), returnData=TRUE)
percentVar <- round(100 * attr(pca, "percentVar"))
p <- ggplot(pca, aes(PC1, PC2, color=cell_line, shape=treatment)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ", percentVar[1], "%")) +
  ylab(paste0("PC2: ", percentVar[2], "%")) +
  theme_bw()
ggsave(filename=file.path(outdir,"figures","PCA_vst.png"), plot=p, width=7, height=5)

# -----------------------------
# Per-cell-line DE
# -----------------------------
cell_lines <- levels(meta$cell_line)
all_res_list <- list()

for (cl in cell_lines) {
  idx <- meta$cell_line == cl
  dds <- DESeqDataSetFromMatrix(
    countData = count_mat[, idx, drop=FALSE],
    colData   = meta[idx, , drop=FALSE],
    design    = ~ treatment
  )
  keep <- rowSums(counts(dds) >= 10) >= 2
  dds <- dds[keep,]
  dds <- DESeq(dds)
  res <- results(dds, contrast=c("treatment","PU001","CTRL"))
  res <- as.data.frame(res)
  res$gene_id <- rownames(res)
  res$gene_name <- gene_name[res$gene_id]
  res <- res[order(res$padj),]
  out_file <- file.path(outdir, paste0("DE_", cl, "_PU001_vs_CTRL.tsv"))
  write.table(res, file=out_file, sep="\t", quote=FALSE, row.names=FALSE)
  all_res_list[[cl]] <- res
}

# -----------------------------
# Hypothesis check: checkpoint genes
# -----------------------------
checkpoint_symbols <- c(
  "HAVCR2",
  "TIGIT",
  "PDCD1",
  "CD80",
  "VSIR",
  "CD86",
  "CD276",
  "LAG3",
  "IDO1",
  "CD274",
  "PDCD1LG2",
  "CTLA4",
  "BTLA"
)

summ <- data.frame()
for (cl in cell_lines) {
  res <- all_res_list[[cl]]
  sub <- res[res$gene_name %in% checkpoint_symbols,
             c("gene_id","gene_name","log2FoldChange","pvalue","padj")]
  sub$cell_line <- cl
  summ <- rbind(summ, sub)
}
glob <- res_all[res_all$gene_name %in% checkpoint_symbols,
                c("gene_id","gene_name","log2FoldChange","pvalue","padj")]
glob$cell_line <- "ALL (cell_line-adjusted)"
summ <- rbind(summ, glob)
summ$supports_hypothesis <- with(summ, !is.na(padj) & padj < 0.05 & log2FoldChange < 0)

write.table(summ, file=file.path(outdir,"checkpoints_summary.tsv"),
            sep="\t", quote=FALSE, row.names=FALSE)

# -----------------------------
# Heatmap
# -----------------------------
heatmap_path <- file.path(outdir, "figures", "checkpoints_heatmap.png")

expr <- assay(vsd)
expr_sym <- gene_name[rownames(expr)]
names(expr_sym) <- rownames(expr)

make_placeholder <- function(msg) {
  png(heatmap_path, width=900, height=650)
  plot.new(); text(0.5, 0.5, msg); dev.off()
}

# 1) find gene_ids for requested symbols
# There can be duplicates — take the first occurrence.
gene_ids_ordered <- sapply(checkpoint_symbols, function(sym) {
  ids <- names(expr_sym)[expr_sym == sym]
  if (length(ids) == 0) return(NA_character_)
  ids[1]
})

present <- !is.na(gene_ids_ordered)
if (sum(present) >= 2) {

  gene_ids_present <- gene_ids_ordered[present]
  symbols_present  <- checkpoint_symbols[present]

  # 2) build matrix in the exact order
  mat <- expr[gene_ids_present, , drop=FALSE]
  rownames(mat) <- symbols_present  # <- sets your desired labels/order

  # 3) drop zero-variance rows
  sds <- apply(mat, 1, sd, na.rm=TRUE)
  keep_rows <- is.finite(sds) & sds > 0
  mat <- mat[keep_rows, , drop=FALSE]

  # 4) align annotation to columns
  ann <- meta[, c("cell_line","treatment"), drop=FALSE]
  ann <- ann[colnames(mat), , drop=FALSE]

  if (nrow(mat) >= 2 && ncol(mat) >= 2) {
    tryCatch({
      pheatmap(
        mat,
        scale = "row",
        cluster_rows = FALSE,
        cluster_cols = TRUE,
        annotation_col = ann,
        show_colnames = FALSE,
        fontsize_row = 10,
        filename = heatmap_path,
        width = 8, height = 6
      )
    }, error = function(e) {
      make_placeholder(paste0("Heatmap failed: ", conditionMessage(e)))
    })
  } else {
    make_placeholder("Not enough variable checkpoint genes for heatmap (after filtering).")
  }

} else {
  make_placeholder("Too few checkpoint genes found for heatmap.")
}

message("Done. Key output: checkpoints_summary.tsv and per-cell-line DE tables.")
