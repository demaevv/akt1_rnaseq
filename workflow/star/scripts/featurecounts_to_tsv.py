#!/usr/bin/env python3
import sys
import os
import pandas as pd

inp = sys.argv[1]
out = sys.argv[2]

df = pd.read_csv(inp, sep="\t", comment="#")

bam_cols = [c for c in df.columns if c.endswith(".bam") or ".bam" in c]
if not bam_cols:
    raise SystemExit("No BAM columns detected in featureCounts output. Check input file format.")

keep = ["Geneid"] + bam_cols
df = df[keep].copy()

def to_sample_id(col: str) -> str:
    b = os.path.basename(col)
    b = b.replace(".Aligned.sortedByCoord.out.bam", "")
    b = b.replace(".Aligned.out.bam", "")
    b = b.replace(".sorted.bam", "")
    b = b.replace(".bam", "")
    return b

new_cols = ["gene_id"] + [to_sample_id(c) for c in bam_cols]
df.columns = new_cols

# sanity
if len(set(df.columns[1:])) != len(df.columns[1:]):
    raise SystemExit("Duplicate sample column names after renaming. Inspect BAM column names.")

df.to_csv(out, sep="\t", index=False)