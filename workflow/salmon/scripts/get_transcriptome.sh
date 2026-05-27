#!/usr/bin/env bash

# Download GENCODE human transcriptome FASTA (and GTF) into ref/
#
# Usage:
#   bash scripts/get_transcriptome.sh              
#   bash scripts/get_transcriptome.sh 49 ref

REL="${1:-49}"
OUTDIR="${2:-ref}"

BASE="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_${REL}"

need() { command -v "$1" >/dev/null 2>&1 || { echo "ERROR: '$1' not found. Install it first."; exit 1; }; }
need aria2c

mkdir -p "$OUTDIR"
cd "$OUTDIR"

DECOMPRESS="gzip -dc"
if command -v pigz >/dev/null 2>&1; then
  DECOMPRESS="pigz -dc"
fi

TX_GZ="gencode.v${REL}.transcripts.fa.gz"
GTF_GZ="gencode.v${REL}.annotation.gtf.gz"

echo "[1/3] Download transcriptome FASTA: $TX_GZ"
aria2c -c -x 8 -s 8 -o "$TX_GZ" "$BASE/$TX_GZ"

echo "[2/3] Download GTF annotation: $GTF_GZ"
aria2c -c -x 8 -s 8 -o "$GTF_GZ" "$BASE/$GTF_GZ"

echo "[3/3] Decompress"
TX_FA="gencode.v${REL}.transcripts.fa"
GTF="gencode.v${REL}.annotation.gtf"

if [[ ! -s "$TX_FA" ]]; then
  $DECOMPRESS "$TX_GZ" > "$TX_FA"
fi
if [[ ! -s "$GTF" ]]; then
  $DECOMPRESS "$GTF_GZ" > "$GTF"
fi

# Stable symlinks for config.yaml
ln -sf "$TX_FA" "gencode.transcriptome.fa"
ln -sf "$GTF"   "gencode.annotation.gtf"

echo "Done:"
echo "  ref/gencode.transcriptome.fa"
echo "  ref/gencode.annotation.gtf"
