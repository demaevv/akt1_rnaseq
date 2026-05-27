#!/usr/bin/env bash

# Usage:
#   bash scripts/get_ref.sh
#   bash scripts/get_ref.sh 49 basic ref

REL="${1:-49}"
ANN="${2:-comprehensive}"
OUTDIR="${3:-ref}"

BASE="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_${REL}"

mkdir -p "$OUTDIR"
cd "$OUTDIR"

# Tools check
need() { command -v "$1" >/dev/null 2>&1 || { echo "ERROR: '$1' not found. Install it first."; exit 1; }; }
need aria2c
need samtools

DECOMPRESS="gzip -dc"
if command -v pigz >/dev/null 2>&1; then
  DECOMPRESS="pigz -dc"
fi

GENOME_GZ="GRCh38.primary_assembly.genome.fa.gz"
if [[ "$ANN" == "basic" ]]; then
  GTF_GZ="gencode.v${REL}.primary_assembly.basic.annotation.gtf.gz"
elif [[ "$ANN" == "comprehensive" ]]; then
  GTF_GZ="gencode.v${REL}.primary_assembly.annotation.gtf.gz"
else
  echo "ERROR: ANN must be 'basic' or 'comprehensive' (got: $ANN)"
  exit 1
fi

echo "[1/4] Download genome FASTA: $GENOME_GZ"
aria2c -c -x 8 -s 8 -o "$GENOME_GZ" "$BASE/$GENOME_GZ"

echo "[2/4] Download GTF annotation: $GTF_GZ"
aria2c -c -x 8 -s 8 -o "$GTF_GZ" "$BASE/$GTF_GZ"

# Optional metadata
for f in "gencode.v${REL}.metadata.HGNC.gz" "gencode.v${REL}.metadata.EntrezGene.gz"; do
  echo "[opt] Download $f (if exists)"
  aria2c -c -x 8 -s 8 -o "$f" "$BASE/$f" || true
done

echo "[3/4] Decompress (creates .fa and .gtf)"
FA="GRCh38.primary_assembly.genome.fa"
GTF="gencode.v${REL}.primary_assembly.${ANN}.annotation.gtf"

if [[ ! -s "$FA" ]]; then
  $DECOMPRESS "$GENOME_GZ" > "$FA"
fi

if [[ ! -s "$GTF" ]]; then
  $DECOMPRESS "$GTF_GZ" > "$GTF"
fi

echo "[4/4] Index FASTA (samtools faidx)"
samtools faidx "$FA"

# Stable symlinks for config.yaml
ln -sf "$FA"  "GRCh38.primary_assembly.genome.fa"
ln -sf "$GTF" "gencode.annotation.gtf"

echo "Done."
echo "ref/GRCh38.primary_assembly.genome.fa"
echo "ref/gencode.annotation.gtf"
