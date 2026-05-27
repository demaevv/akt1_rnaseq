#!/usr/bin/env python3
"""
Pathway/Enrichment Analysis for THP-1: PU-001 vs Control
Master's thesis — ITMO University
Uses STAR-based DESeq2 output.
"""

import os
import warnings
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import seaborn as sns
from scipy import stats
from scipy.stats import hypergeom
from statsmodels.stats.multitest import multipletests
import gseapy as gp
from gseapy import gseaplot

warnings.filterwarnings('ignore')
sns.set_style('whitegrid')
sns.set_context('paper', font_scale=1.1)

# ============================================================
# CONFIGURATION
# ============================================================
INPUT_TSV = 'DE_THP-1_PU001_vs_CTRL.tsv'
OUTDIR = 'pathway_analysis_THP1'
PADJ_CUTOFF = 0.05
LFC_CUTOFF = 1.0
TOP_N = 200
GSEA_PERMUTATIONS = 1000
GSEA_MIN_SIZE = 15
GSEA_MAX_SIZE = 500
SEED = 42

os.makedirs(OUTDIR, exist_ok=True)
os.makedirs(f'{OUTDIR}/figures', exist_ok=True)
os.makedirs(f'{OUTDIR}/tables', exist_ok=True)
os.makedirs(f'{OUTDIR}/gsea_output', exist_ok=True)

# ============================================================
# STEP 1: DATA LOADING AND PREPARATION
# ============================================================
print("=" * 60)
print("STEP 1: Data Loading and Filtering")
print("=" * 60)

df = pd.read_csv(INPUT_TSV, sep='\t')
print(f"Total genes in DE table: {len(df)}")
print(f"Columns: {list(df.columns)}")

# Remove genes with NA padj or pvalue
df_clean = df.dropna(subset=['padj', 'pvalue', 'log2FoldChange', 'gene_name']).copy()
print(f"Genes with non-NA padj: {len(df_clean)}")

# Filter: padj < 0.05 AND |log2FC| >= 1
sig = df_clean[(df_clean['padj'] < PADJ_CUTOFF) & 
               (df_clean['log2FoldChange'].abs() >= LFC_CUTOFF)].copy()
print(f"\nGenes passing filter (padj<{PADJ_CUTOFF}, |log2FC|>={LFC_CUTOFF}): {len(sig)}")

# Split UP and DOWN
sig_up = sig[sig['log2FoldChange'] >= LFC_CUTOFF].copy()
sig_down = sig[sig['log2FoldChange'] <= -LFC_CUTOFF].copy()
print(f"  UP-regulated genes: {len(sig_up)}")
print(f"  DOWN-regulated genes: {len(sig_down)}")

# Rank by ascending padj and take top 200
sig_up = sig_up.sort_values('padj', ascending=True)
sig_down = sig_down.sort_values('padj', ascending=True)

n_up_before = len(sig_up)
n_down_before = len(sig_down)

top_up = sig_up.head(TOP_N).copy()
top_down = sig_down.head(TOP_N).copy()

print(f"\n  UP genes before top-{TOP_N} cutoff: {n_up_before}" +
      (f" (using all {n_up_before})" if n_up_before < TOP_N else f" (taking top {TOP_N})"))
print(f"  DOWN genes before top-{TOP_N} cutoff: {n_down_before}" +
      (f" (using all {n_down_before})" if n_down_before < TOP_N else f" (taking top {TOP_N})"))

up_genes = top_up['gene_name'].str.upper().tolist()
down_genes = top_down['gene_name'].str.upper().tolist()

# Save gene lists
top_up.to_csv(f'{OUTDIR}/tables/top200_UP_genes.tsv', sep='\t', index=False)
top_down.to_csv(f'{OUTDIR}/tables/top200_DOWN_genes.tsv', sep='\t', index=False)

# Background gene set (all expressed genes with non-NA values)
background_genes = df_clean['gene_name'].str.upper().dropna().unique().tolist()
print(f"  Background gene universe size: {len(background_genes)}")

# Summary statistics
print(f"\n  UP genes — padj range: [{top_up['padj'].min():.2e}, {top_up['padj'].max():.2e}]")
print(f"  UP genes — log2FC range: [{top_up['log2FoldChange'].min():.2f}, {top_up['log2FoldChange'].max():.2f}]")
print(f"  DOWN genes — padj range: [{top_down['padj'].min():.2e}, {top_down['padj'].max():.2e}]")
print(f"  DOWN genes — log2FC range: [{top_down['log2FoldChange'].min():.2f}, {top_down['log2FoldChange'].max():.2f}]")

# ============================================================
# STEP 2A: MSigDB-based ORA (Hypergeometric Test)
# ============================================================
print("\n" + "=" * 60)
print("STEP 2A: MSigDB Over-Representation Analysis (ORA)")
print("=" * 60)

def run_ora_msigdb(gene_list, gene_set_name, background, label=""):
    """Run ORA using gseapy's enrich function against a gene set library."""
    try:
        enr = gp.enrich(
            gene_list=gene_list,
            gene_sets=gene_set_name,
            background=background,
            outdir=None,
            no_plot=True,
            cutoff=1.0,  # get all results, filter later
            verbose=False
        )
        res = enr.res2d.copy()
        if len(res) > 0:
            res = res.sort_values('Adjusted P-value', ascending=True)
        return res
    except Exception as e:
        print(f"  Warning: ORA failed for {gene_set_name} ({label}): {e}")
        return pd.DataFrame()

# MSigDB collections to query
# Note: Download .gmt files from https://www.gsea-msigdb.org/gsea/msigdb/
# Or use gseapy Enrichr library names as approximation
# For precise MSigDB analysis, use downloaded .gmt files

# Try to use Msigdb class first, fall back to Enrichr libraries
msigdb_collections = {}

# Try downloading GMT files via gseapy's Msigdb class
try:
    from gseapy import Msigdb
    msig = Msigdb()
    # Try getting gene sets
    # If this doesn't work, fall back to Enrichr-based libraries
    raise NotImplementedError("Falling back to GMT file approach")
except:
    # Use Enrichr library names or download GMT files
    # Check if GMT files exist locally; if not, try to download
    GMT_DIR = f'{OUTDIR}/gmt_files'
    os.makedirs(GMT_DIR, exist_ok=True)
    
    # First, try using gseapy's get_library for Enrichr-hosted MSigDB sets
    print("Using gseapy Enrichr-hosted libraries for ORA...")
    
    msigdb_map = {
        'H_Hallmark': 'MSigDB_Hallmark_2020',
        'C6_Oncogenic': 'MSigDB_Oncogenic_Signatures',
    }
    
    # For C2:CP and C7, we need .gmt files or alternative approaches
    # Try to download from MSigDB directly
    import urllib.request
    
    gmt_urls = {
        'H_Hallmark': 'https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2024.1.Hs/h.all.v2024.1.Hs.symbols.gmt',
        'C2_CP': 'https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2024.1.Hs/c2.cp.v2024.1.Hs.symbols.gmt',
        'C6_Oncogenic': 'https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2024.1.Hs/c6.all.v2024.1.Hs.symbols.gmt',
        'C7_Immunologic': 'https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2024.1.Hs/c7.all.v2024.1.Hs.symbols.gmt',
    }
    
    for name, url in gmt_urls.items():
        gmt_path = f'{GMT_DIR}/{name}.gmt'
        if not os.path.exists(gmt_path):
            try:
                print(f"  Downloading {name} GMT file...")
                urllib.request.urlretrieve(url, gmt_path)
                msigdb_collections[name] = gmt_path
            except Exception as e:
                print(f"  Could not download {name}: {e}")
                # Fall back to Enrichr library
                if name in msigdb_map:
                    msigdb_collections[name] = msigdb_map[name]
        else:
            msigdb_collections[name] = gmt_path

    # If GMT download failed, use Enrichr libraries as fallback
    if 'H_Hallmark' not in msigdb_collections:
        msigdb_collections['H_Hallmark'] = 'MSigDB_Hallmark_2020'
    if 'C6_Oncogenic' not in msigdb_collections:
        msigdb_collections['C6_Oncogenic'] = 'MSigDB_Oncogenic_Signatures'
    # For C2:CP, try individual Enrichr libraries
    if 'C2_CP' not in msigdb_collections:
        msigdb_collections['C2_CP_Reactome'] = 'Reactome_2022'
        msigdb_collections['C2_CP_KEGG'] = 'KEGG_2021_Human'
    if 'C7_Immunologic' not in msigdb_collections:
        # C7 is not directly available via Enrichr; try alternate
        print("  Note: C7 Immunologic signatures not available via Enrichr; will skip or use GMT")

# Run ORA for each collection, for both UP and DOWN gene lists
ora_results = {}
for direction, gene_list in [('UP', up_genes), ('DOWN', down_genes)]:
    print(f"\n--- ORA for {direction}-regulated genes ({len(gene_list)} genes) ---")
    for coll_name, coll_source in msigdb_collections.items():
        print(f"  Running ORA: {coll_name}...")
        res = run_ora_msigdb(gene_list, coll_source, background_genes, f"{direction}_{coll_name}")
        ora_results[f'{direction}_{coll_name}'] = res
        if len(res) > 0:
            top10 = res.head(10)
            print(f"    Top pathway: {top10.iloc[0]['Term'] if len(top10) > 0 else 'N/A'}")
            print(f"    Significant terms (Adj.P<0.05): {(res['Adjusted P-value'] < 0.05).sum()}")
        # Save full results
        res.to_csv(f'{OUTDIR}/tables/ORA_{direction}_{coll_name}.tsv', sep='\t', index=False)

# ============================================================
# STEP 2B: Enrichr-based Enrichment (Cross-validation)
# ============================================================
print("\n" + "=" * 60)
print("STEP 2B: Enrichr-based Enrichment Analysis")
print("=" * 60)

enrichr_libraries = [
    'Reactome_2022',
    'KEGG_2021_Human',
    'GO_Biological_Process_2023',
    'MSigDB_Hallmark_2020',
    'ChEA_2022',
    'TRRUST_Transcription_Factors_2019'
]

enrichr_results = {}
for direction, gene_list in [('UP', up_genes), ('DOWN', down_genes)]:
    print(f"\n--- Enrichr for {direction}-regulated genes ---")
    for lib in enrichr_libraries:
        print(f"  Querying Enrichr: {lib}...")
        try:
            enr = gp.enrichr(
                gene_list=gene_list,
                gene_sets=lib,
                organism='hsapiens',
                outdir=None,
                no_plot=True,
                cutoff=1.0,
                verbose=False
            )
            res = enr.res2d.copy()
            enrichr_results[f'{direction}_{lib}'] = res
            if len(res) > 0:
                sig_count = (res['Adjusted P-value'] < 0.05).sum()
                print(f"    Significant terms (Adj.P<0.05): {sig_count}")
                if sig_count > 0:
                    print(f"    Top term: {res.iloc[0]['Term']}")
            res.to_csv(f'{OUTDIR}/tables/Enrichr_{direction}_{lib}.tsv', sep='\t', index=False)
        except Exception as e:
            print(f"    Error: {e}")
            enrichr_results[f'{direction}_{lib}'] = pd.DataFrame()

# ============================================================
# STEP 2C: Pre-ranked GSEA (Full Gene List)
# ============================================================
print("\n" + "=" * 60)
print("STEP 2C: Pre-ranked GSEA")
print("=" * 60)

# Prepare ranking metric: sign(log2FC) * (-log10(pvalue))
df_gsea = df_clean[['gene_name', 'log2FoldChange', 'pvalue']].dropna().copy()
df_gsea['gene_name'] = df_gsea['gene_name'].str.upper()

# Handle p-value = 0 by clipping to minimum non-zero value
min_pval = df_gsea.loc[df_gsea['pvalue'] > 0, 'pvalue'].min()
df_gsea['pvalue_clipped'] = df_gsea['pvalue'].clip(lower=min_pval * 0.1)
df_gsea['rank_metric'] = np.sign(df_gsea['log2FoldChange']) * (-np.log10(df_gsea['pvalue_clipped']))

# Remove duplicates (keep highest absolute rank metric)
df_gsea['abs_rank'] = df_gsea['rank_metric'].abs()
df_gsea = df_gsea.sort_values('abs_rank', ascending=False).drop_duplicates(subset='gene_name', keep='first')
df_gsea = df_gsea.sort_values('rank_metric', ascending=False).reset_index(drop=True)

# Prepare rnk dataframe
rnk = df_gsea[['gene_name', 'rank_metric']].copy()
rnk.columns = ['gene_name', 'rank_metric']
print(f"Genes in ranked list: {len(rnk)}")
print(f"Rank metric range: [{rnk['rank_metric'].min():.2f}, {rnk['rank_metric'].max():.2f}]")

# Save ranked list
rnk.to_csv(f'{OUTDIR}/tables/ranked_gene_list.rnk', sep='\t', index=False, header=False)

# Run GSEA against each collection
gsea_results = {}
gsea_objects = {}

for coll_name, coll_source in msigdb_collections.items():
    print(f"\n  Running GSEA: {coll_name}...")
    try:
        pre_res = gp.prerank(
            rnk=rnk,
            gene_sets=coll_source,
            permutation_num=GSEA_PERMUTATIONS,
            min_size=GSEA_MIN_SIZE,
            max_size=GSEA_MAX_SIZE,
            weight=1.0,
            ascending=False,
            threads=4,
            outdir=f'{OUTDIR}/gsea_output/{coll_name}',
            format='png',
            seed=SEED,
            verbose=False
        )
        res = pre_res.res2d.copy()
        gsea_results[coll_name] = res
        gsea_objects[coll_name] = pre_res
        
        # Report top positively and negatively enriched
        pos = res[res['NES'] > 0].sort_values('FDR q-val', ascending=True)
        neg = res[res['NES'] < 0].sort_values('FDR q-val', ascending=True)
        print(f"    Positively enriched (FDR<0.25): {(pos['FDR q-val'] < 0.25).sum()}")
        print(f"    Negatively enriched (FDR<0.25): {(neg['FDR q-val'] < 0.25).sum()}")
        if len(pos) > 0:
            print(f"    Top positive: {pos.iloc[0]['Term']} (NES={pos.iloc[0]['NES']:.2f}, FDR={pos.iloc[0]['FDR q-val']:.4f})")
        if len(neg) > 0:
            print(f"    Top negative: {neg.iloc[0]['Term']} (NES={neg.iloc[0]['NES']:.2f}, FDR={neg.iloc[0]['FDR q-val']:.4f})")
        
        res.to_csv(f'{OUTDIR}/tables/GSEA_{coll_name}.tsv', sep='\t', index=False)
    except Exception as e:
        print(f"    Error: {e}")

# ============================================================
# STEP 3: VISUALIZATION
# ============================================================
print("\n" + "=" * 60)
print("STEP 3: Generating Figures")
print("=" * 60)

def get_best_gsea_for_direction(gsea_results_dict, direction='positive'):
    """Find the single most significant pathway across all GSEA collections for a direction."""
    best_term = None
    best_fdr = 1.0
    best_overlap = 0
    best_collection = None
    
    for coll_name, res in gsea_results_dict.items():
        if direction == 'positive':
            subset = res[res['NES'] > 0].copy()
        else:
            subset = res[res['NES'] < 0].copy()
        
        if len(subset) == 0:
            continue
            
        subset = subset.sort_values(['FDR q-val', 'NES'], ascending=[True, direction != 'positive'])
        top = subset.iloc[0]
        
        # Parse leading edge size
        try:
            le_size = len(str(top.get('Lead_genes', top.get('lead_genes', ''))).split(';'))
        except:
            le_size = 0
        
        if top['FDR q-val'] < best_fdr or (top['FDR q-val'] == best_fdr and le_size > best_overlap):
            best_fdr = top['FDR q-val']
            best_term = top['Term']
            best_overlap = le_size
            best_collection = coll_name
    
    return best_term, best_collection, best_fdr, best_overlap

def make_gsea_running_plot(gsea_obj, term_name, title_prefix, filename):
    """Create a classic GSEA running enrichment score plot."""
    try:
        # Use gseapy's built-in plot
        axs = gsea_obj.plot(terms=[term_name], show_ranking=True, figsize=(8, 6))
        if axs is not None:
            fig = axs[0].figure if hasattr(axs, '__iter__') else axs.figure
            fig.savefig(filename, dpi=300, bbox_inches='tight')
            plt.close(fig)
            print(f"  Saved: {filename}")
            return
    except Exception as e:
        print(f"  Built-in plot failed ({e}), trying manual approach...")
    
    # Manual GSEA plot
    try:
        result = gsea_obj.results[term_name]
        ranking = gsea_obj.ranking
        
        fig = plt.figure(figsize=(8, 6))
        gs = gridspec.GridSpec(3, 1, height_ratios=[3, 0.5, 1.5], hspace=0.05)
        
        # Panel 1: Running ES
        ax1 = fig.add_subplot(gs[0])
        es_values = result.get('RES', result.get('es', []))
        if hasattr(es_values, '__len__') and len(es_values) > 0:
            ax1.plot(range(len(es_values)), es_values, color='green', linewidth=2)
            ax1.axhline(y=0, color='grey', linestyle='--', linewidth=0.5)
            ax1.set_ylabel('Enrichment Score (ES)')
            
            # Get NES and FDR from res2d
            res_row = gsea_obj.res2d[gsea_obj.res2d['Term'] == term_name].iloc[0]
            nes = res_row['NES']
            fdr = res_row['FDR q-val']
            
            ax1.set_title(f'{title_prefix}\n{term_name}\nNES={nes:.2f}, FDR={fdr:.4f}', fontsize=10)
        
        # Panel 2: Hit positions
        ax2 = fig.add_subplot(gs[1], sharex=ax1)
        hits = result.get('hits', [])
        if hasattr(hits, '__len__') and len(hits) > 0:
            for hit in hits:
                ax2.axvline(x=hit, color='black', linewidth=0.3)
        ax2.set_ylabel('')
        ax2.set_yticks([])
        
        # Panel 3: Ranked metric
        ax3 = fig.add_subplot(gs[2], sharex=ax1)
        if hasattr(ranking, '__len__') and len(ranking) > 0:
            ax3.fill_between(range(len(ranking)), ranking, color='lightgrey')
            ax3.axhline(y=0, color='grey', linestyle='--', linewidth=0.5)
            ax3.set_ylabel('Ranked Metric')
            ax3.set_xlabel('Gene Rank')
        
        plt.tight_layout()
        fig.savefig(filename, dpi=300, bbox_inches='tight')
        plt.close(fig)
        print(f"  Saved: {filename}")
    except Exception as e:
        print(f"  Manual GSEA plot also failed: {e}")

def make_dotplot(results_df, title, filename, top_n=10, metric_col='NES', 
                 size_col=None, color_col='FDR q-val'):
    """Create a dotplot for enrichment results."""
    fig, ax = plt.subplots(figsize=(10, 6))
    
    if len(results_df) == 0:
        ax.text(0.5, 0.5, 'No significant results', ha='center', va='center', fontsize=14)
        fig.savefig(filename, dpi=300, bbox_inches='tight')
        plt.close(fig)
        return
    
    data = results_df.head(top_n).copy()
    data = data.iloc[::-1]  # Reverse for bottom-to-top display
    
    # Determine if this is GSEA or ORA/Enrichr results
    if 'NES' in data.columns:
        x_values = data['NES'].astype(float)
        x_label = 'Normalized Enrichment Score (NES)'
        # Parse overlap count from lead_genes
        if 'Lead_genes' in data.columns:
            sizes = data['Lead_genes'].apply(lambda x: len(str(x).split(';')) if pd.notna(x) else 5)
        elif 'lead_genes' in data.columns:
            sizes = data['lead_genes'].apply(lambda x: len(str(x).split(';')) if pd.notna(x) else 5)
        else:
            sizes = [30] * len(data)
        fdr_values = -np.log10(data[color_col].astype(float).clip(lower=1e-50))
    elif 'Combined Score' in data.columns:
        x_values = data['Combined Score'].astype(float)
        x_label = 'Combined Score'
        if 'Overlap' in data.columns:
            sizes = data['Overlap'].apply(lambda x: int(str(x).split('/')[0]) if '/' in str(x) else 5)
        else:
            sizes = [30] * len(data)
        fdr_values = -np.log10(data['Adjusted P-value'].astype(float).clip(lower=1e-50))
    else:
        x_values = -np.log10(data['Adjusted P-value'].astype(float).clip(lower=1e-50))
        x_label = '-log10(Adjusted P-value)'
        sizes = [30] * len(data)
        fdr_values = x_values
    
    # Truncate long pathway names
    y_labels = [name[:60] + '...' if len(str(name)) > 60 else str(name) for name in data['Term']]
    
    scatter = ax.scatter(x_values, range(len(data)), 
                         c=fdr_values, cmap='viridis_r',
                         s=np.array(sizes) * 8 + 20,
                         edgecolors='black', linewidth=0.5, alpha=0.8)
    
    ax.set_yticks(range(len(data)))
    ax.set_yticklabels(y_labels, fontsize=8)
    ax.set_xlabel(x_label)
    ax.set_title(title, fontsize=11, fontweight='bold')
    
    plt.colorbar(scatter, ax=ax, label='-log10(FDR)', shrink=0.8)
    
    # Add size legend
    legend_sizes = [5, 20, 50]
    for s in legend_sizes:
        ax.scatter([], [], s=s*8+20, c='grey', edgecolors='black', linewidth=0.5,
                   label=f'{s} genes')
    ax.legend(title='Overlap', loc='lower right', fontsize=7)
    
    plt.tight_layout()
    fig.savefig(filename, dpi=300, bbox_inches='tight')
    plt.close(fig)
    print(f"  Saved: {filename}")

def make_barplot(results_df, title, filename, top_n=10):
    """Create a barplot for enrichment results."""
    fig, ax = plt.subplots(figsize=(10, 6))
    
    if len(results_df) == 0:
        ax.text(0.5, 0.5, 'No significant results', ha='center', va='center', fontsize=14)
        fig.savefig(filename, dpi=300, bbox_inches='tight')
        plt.close(fig)
        return
    
    data = results_df.head(top_n).copy()
    data = data.iloc[::-1]
    
    if 'NES' in data.columns:
        x_values = data['NES'].astype(float)
        x_label = 'Normalized Enrichment Score (NES)'
        colors = ['#E74C3C' if v > 0 else '#3498DB' for v in x_values]
    elif 'Adjusted P-value' in data.columns:
        x_values = -np.log10(data['Adjusted P-value'].astype(float).clip(lower=1e-50))
        x_label = '-log10(Adjusted P-value)'
        colors = '#E74C3C'
    else:
        x_values = range(len(data))
        x_label = 'Score'
        colors = '#E74C3C'
    
    y_labels = [name[:60] + '...' if len(str(name)) > 60 else str(name) for name in data['Term']]
    
    ax.barh(range(len(data)), x_values, color=colors, edgecolor='black', linewidth=0.5, alpha=0.85)
    ax.set_yticks(range(len(data)))
    ax.set_yticklabels(y_labels, fontsize=8)
    ax.set_xlabel(x_label)
    ax.set_title(title, fontsize=11, fontweight='bold')
    
    if 'NES' in data.columns:
        ax.axvline(x=0, color='black', linewidth=0.8)
    
    plt.tight_layout()
    fig.savefig(filename, dpi=300, bbox_inches='tight')
    plt.close(fig)
    print(f"  Saved: {filename}")

# --- Generate UP-regulated figures ---
print("\nGenerating UP-regulated figures...")

# Find best GSEA pathway for UP (positive NES)
best_up_term, best_up_coll, best_up_fdr, best_up_overlap = get_best_gsea_for_direction(gsea_results, 'positive')

if best_up_term and best_up_coll in gsea_objects:
    # Figure 1: GSEA running plot for most significant UP pathway
    make_gsea_running_plot(
        gsea_objects[best_up_coll], best_up_term,
        "UP-regulated",
        f'{OUTDIR}/figures/Fig1_GSEA_running_UP.png'
    )

# For dotplot and barplot, use Hallmark GSEA results (positive NES)
hallmark_key = 'H_Hallmark'
if hallmark_key in gsea_results and len(gsea_results[hallmark_key]) > 0:
    up_gsea = gsea_results[hallmark_key][gsea_results[hallmark_key]['NES'] > 0].sort_values('FDR q-val')
    make_dotplot(up_gsea, 'UP-regulated: Top Hallmark Pathways (GSEA)', 
                 f'{OUTDIR}/figures/Fig2_dotplot_UP.png')
    make_barplot(up_gsea, 'UP-regulated: Top Hallmark Pathways (GSEA)', 
                 f'{OUTDIR}/figures/Fig3_barplot_UP.png')
else:
    # Fall back to ORA/Enrichr results
    for key in ora_results:
        if key.startswith('UP_') and len(ora_results[key]) > 0:
            make_dotplot(ora_results[key], f'UP-regulated: {key} (ORA)',
                         f'{OUTDIR}/figures/Fig2_dotplot_UP.png')
            make_barplot(ora_results[key], f'UP-regulated: {key} (ORA)',
                         f'{OUTDIR}/figures/Fig3_barplot_UP.png')
            break

# --- Generate DOWN-regulated figures ---
print("\nGenerating DOWN-regulated figures...")

best_dn_term, best_dn_coll, best_dn_fdr, best_dn_overlap = get_best_gsea_for_direction(gsea_results, 'negative')

if best_dn_term and best_dn_coll in gsea_objects:
    make_gsea_running_plot(
        gsea_objects[best_dn_coll], best_dn_term,
        "DOWN-regulated",
        f'{OUTDIR}/figures/Fig4_GSEA_running_DOWN.png'
    )

if hallmark_key in gsea_results and len(gsea_results[hallmark_key]) > 0:
    dn_gsea = gsea_results[hallmark_key][gsea_results[hallmark_key]['NES'] < 0].sort_values('FDR q-val')
    make_dotplot(dn_gsea, 'DOWN-regulated: Top Hallmark Pathways (GSEA)',
                 f'{OUTDIR}/figures/Fig5_dotplot_DOWN.png')
    make_barplot(dn_gsea, 'DOWN-regulated: Top Hallmark Pathways (GSEA)',
                 f'{OUTDIR}/figures/Fig6_barplot_DOWN.png')
else:
    for key in ora_results:
        if key.startswith('DOWN_') and len(ora_results[key]) > 0:
            make_dotplot(ora_results[key], f'DOWN-regulated: {key} (ORA)',
                         f'{OUTDIR}/figures/Fig5_dotplot_DOWN.png')
            make_barplot(ora_results[key], f'DOWN-regulated: {key} (ORA)',
                         f'{OUTDIR}/figures/Fig6_barplot_DOWN.png')
            break

# ============================================================
# STEP 5: GENERATE MARKDOWN REPORT
# ============================================================
print("\n" + "=" * 60)
print("STEP 5: Generating Markdown Report")
print("=" * 60)

def format_ora_table(res_df, top_n=10):
    """Format ORA results as markdown table."""
    if len(res_df) == 0:
        return "No significant results found.\n"
    data = res_df.head(top_n)
    lines = ["| Rank | Pathway | Overlap | P-value | FDR q-value | Genes |",
             "|------|---------|---------|---------|-------------|-------|"]
    for i, (_, row) in enumerate(data.iterrows(), 1):
        term = str(row.get('Term', 'N/A'))[:80]
        overlap = str(row.get('Overlap', 'N/A'))
        pval = f"{float(row.get('P-value', 1)):.2e}"
        adjp = f"{float(row.get('Adjusted P-value', 1)):.2e}"
        genes = str(row.get('Genes', ''))[:100]
        lines.append(f"| {i} | {term} | {overlap} | {pval} | {adjp} | {genes} |")
    return '\n'.join(lines) + '\n'

def format_enrichr_table(res_df, top_n=10):
    """Format Enrichr results as markdown table."""
    if len(res_df) == 0:
        return "No significant results found.\n"
    data = res_df.head(top_n)
    lines = ["| Rank | Term | Overlap | Adj. P-value | Combined Score | Genes |",
             "|------|------|---------|--------------|----------------|-------|"]
    for i, (_, row) in enumerate(data.iterrows(), 1):
        term = str(row.get('Term', 'N/A'))[:80]
        overlap = str(row.get('Overlap', 'N/A'))
        adjp = f"{float(row.get('Adjusted P-value', 1)):.2e}"
        cscore = f"{float(row.get('Combined Score', 0)):.1f}"
        genes = str(row.get('Genes', ''))[:100]
        lines.append(f"| {i} | {term} | {overlap} | {adjp} | {cscore} | {genes} |")
    return '\n'.join(lines) + '\n'

def format_gsea_table(res_df, direction='positive', top_n=10):
    """Format GSEA results as markdown table."""
    if len(res_df) == 0:
        return "No significant results found.\n"
    if direction == 'positive':
        data = res_df[res_df['NES'] > 0].sort_values('FDR q-val').head(top_n)
    else:
        data = res_df[res_df['NES'] < 0].sort_values('FDR q-val').head(top_n)
    
    if len(data) == 0:
        return "No enriched pathways in this direction.\n"
    
    lines = ["| Rank | Pathway | NES | FDR q-val | Leading Edge Size |",
             "|------|---------|-----|-----------|-------------------|"]
    for i, (_, row) in enumerate(data.iterrows(), 1):
        term = str(row.get('Term', 'N/A'))[:80]
        nes = f"{float(row['NES']):.2f}"
        fdr = f"{float(row['FDR q-val']):.4f}"
        # Parse leading edge
        le_genes = str(row.get('Lead_genes', row.get('lead_genes', '')))
        le_size = len(le_genes.split(';')) if le_genes and le_genes != 'nan' else 0
        lines.append(f"| {i} | {term} | {nes} | {fdr} | {le_size} |")
    return '\n'.join(lines) + '\n'

# Build the Markdown report
report = f"""# Pathway and Enrichment Analysis Report
## THP-1 Cell Line: PU-001 (200 μM, 24h) vs Control
### STAR + featureCounts → DESeq2 Pipeline

---

## 1. Data Filtering Statistics

- **Total genes in DE table:** {len(df)}
- **Genes with non-NA padj:** {len(df_clean)}
- **Genes passing filter (padj<{PADJ_CUTOFF}, |log2FC|≥{LFC_CUTOFF}):** {len(sig)}
  - **UP-regulated:** {n_up_before}{' (all used)' if n_up_before < TOP_N else f' (top {TOP_N} used)'}
  - **DOWN-regulated:** {n_down_before}{' (all used)' if n_down_before < TOP_N else f' (top {TOP_N} used)'}
- **Background universe:** {len(background_genes)} genes

### Top 10 UP-regulated genes (by padj):
{top_up[['gene_name','log2FoldChange','padj']].head(10).to_markdown(index=False)}

### Top 10 DOWN-regulated genes (by padj):
{top_down[['gene_name','log2FoldChange','padj']].head(10).to_markdown(index=False)}

---

## 2A. MSigDB Over-Representation Analysis (ORA)

"""

for direction in ['UP', 'DOWN']:
    report += f"### {direction}-regulated genes\n\n"
    for coll_name in msigdb_collections:
        key = f'{direction}_{coll_name}'
        if key in ora_results:
            report += f"#### {coll_name}\n\n"
            report += format_ora_table(ora_results[key])
            report += "\n"

report += """
---

## 2B. Enrichr-based Enrichment (Cross-validation)

"""

for direction in ['UP', 'DOWN']:
    report += f"### {direction}-regulated genes\n\n"
    for lib in enrichr_libraries:
        key = f'{direction}_{lib}'
        if key in enrichr_results:
            report += f"#### {lib}\n\n"
            report += format_enrichr_table(enrichr_results[key])
            report += "\n"

report += """
---

## 2C. Pre-ranked GSEA (Full Gene List)

"""

for coll_name in msigdb_collections:
    if coll_name in gsea_results:
        report += f"### {coll_name}\n\n"
        report += "#### Top 10 Positively Enriched Pathways\n\n"
        report += format_gsea_table(gsea_results[coll_name], 'positive')
        report += "\n#### Top 10 Negatively Enriched Pathways\n\n"
        report += format_gsea_table(gsea_results[coll_name], 'negative')
        report += "\n"

report += f"""
---

## 3. Figures

### UP-regulated analysis

![GSEA Running Plot — Most significant UP pathway](figures/Fig1_GSEA_running_UP.png)

![Dotplot — Top 10 UP-enriched pathways](figures/Fig2_dotplot_UP.png)

![Barplot — Top 10 UP-enriched pathways](figures/Fig3_barplot_UP.png)

### DOWN-regulated analysis

![GSEA Running Plot — Most significant DOWN pathway](figures/Fig4_GSEA_running_DOWN.png)

![Dotplot — Top 10 DOWN-enriched pathways](figures/Fig5_dotplot_DOWN.png)

![Barplot — Top 10 DOWN-enriched pathways](figures/Fig6_barplot_DOWN.png)

---

## 4. Biological Interpretation

*See the companion interpretive discussion document.*

---

*Generated by pathway_analysis_thp1.py*
*Analysis date: {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M')}*
"""

# Write report
report_path = f'{OUTDIR}/pathway_analysis_report.md'
with open(report_path, 'w') as f:
    f.write(report)

print(f"\nReport saved to: {report_path}")
print(f"Figures saved in: {OUTDIR}/figures/")
print(f"Tables saved in: {OUTDIR}/tables/")
print("\nDone!")