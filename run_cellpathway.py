#!/usr/bin/env python
"""Example: run the full CellPathway pipeline on the bundled autism dataset.

Usage
-----
    cd CellPathway
    python run_cellpathway.py

This script mirrors the notebook workflow (Steps 1-8) using the Python API.
"""

import os
import pandas as pd
from cellpathway import CellPathway
from cellpathway.tad import annotate_tad

# ------------------------------------------------------------------
# 1. Paths — adjust if your data lives elsewhere
# ------------------------------------------------------------------
BASE_DIR = os.path.dirname(os.path.abspath(__file__))

ENHANCER_DIR = os.path.join(BASE_DIR, "data", "Atlas")
DNM_FILE = os.path.join(BASE_DIR, "example", "autism_dnm.txt")
OUTPUT_DIR = os.path.join(BASE_DIR, "example")
TAD_FILE = os.path.join(BASE_DIR, "data", "tad_w_boundary_08.bed")
ELEMENTS_BB = os.path.join(BASE_DIR, "data", "genes_w_noncoding.bb")
SFARI_FILE = os.path.join(BASE_DIR, "example", "SFARI_Gene.csv")

CADD_THRESHOLD = 10

# ------------------------------------------------------------------
# 2. Run enrichment analysis (Steps 1-7)
# ------------------------------------------------------------------
cp = CellPathway(
    enhancer_dir=ENHANCER_DIR,
    dnm_file=DNM_FILE,
    output_dir=OUTPUT_DIR,
    cadd_threshold=CADD_THRESHOLD,
)
results = cp.run()

# Show top 10 enriched cell types
top10 = results.sort_values("P_FDR").head(10)
print("\n=== Top 10 enriched cell types ===")
print(top10[["cell", "Fold_enrichment", "P_FDR"]].to_string(index=False))

# ------------------------------------------------------------------
# 3. TAD annotation (Step 8) — for the top cell type
# ------------------------------------------------------------------
top_cell = top10.iloc[0]["cell"]
print(f"\nRunning TAD annotation for: {top_cell}")

overlap_bed = os.path.join(
    OUTPUT_DIR,
    f"dnm_enhc_overlap_cadd_{CADD_THRESHOLD}",
    f"{top_cell}_dnm.bed",
)

# Load SFARI gene list
sfari = pd.read_csv(SFARI_FILE)
sfari_genes = sfari["gene-symbol"].tolist()

tad_results = annotate_tad(
    overlap_bed_path=overlap_bed,
    tad_path=TAD_FILE,
    elements_bb_path=ELEMENTS_BB,
    gene_list=sfari_genes,
    output_path=os.path.join(OUTPUT_DIR, f"tad_{top_cell}_autism.csv"),
)

print(f"\nTAD annotation rows: {len(tad_results)}")
print("Sample output:")
print(tad_results[["enh_chr", "enh_start", "enh_end", "n_genes",
                    "known_gene_overlap"]].head(10).to_string(index=False))
