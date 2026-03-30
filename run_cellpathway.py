#!/usr/bin/env python
"""CellPathway: cell type-specific enhancer enrichment analysis.

Usage
-----
    # Step 1-7: Enrichment analysis
    python run_cellpathway.py enrich \
        --enhancer-dir data/Atlas \
        --dnm-file example/autism_dnm.txt \
        --output-dir example \
        --cadd-threshold 10

    # Step 8: TAD annotation (run after enrichment)
    python run_cellpathway.py tad \
        --overlap-bed example/dnm_enhc_overlap_cadd_10/Fetal_brain_dnm.bed \
        --tad-file data/tad_w_boundary_08.bed \
        --elements-bb data/genes_w_noncoding.bb \
        --gene-list example/SFARI_Gene.csv \
        --output example/tad_Fetal_brain_autism.csv
"""

import sys
import os
import argparse

# Allow importing cellpathway without pip install
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from cellpathway.core import CellPathway
from cellpathway.tad import annotate_tad


def cmd_enrich(args):
    cp = CellPathway(
        enhancer_dir=args.enhancer_dir,
        dnm_file=args.dnm_file,
        output_dir=args.output_dir,
        cadd_threshold=args.cadd_threshold,
        bedtools_path=args.bedtools_path,
    )
    results = cp.run()

    top = results.sort_values("P_FDR").head(10)
    print("\nTop 10 enriched cell types:")
    print(top[["cell", "Fold_enrichment", "P_FDR"]].to_string(index=False))


def cmd_tad(args):
    import pandas as pd

    gene_list = None
    if args.gene_list:
        genes_df = pd.read_csv(args.gene_list)
        gene_list = genes_df["gene-symbol"].tolist()

    results = annotate_tad(
        overlap_bed_path=args.overlap_bed,
        tad_path=args.tad_file,
        elements_bb_path=args.elements_bb,
        gene_list=gene_list,
        output_path=args.output,
    )
    print(f"\nTAD annotation rows: {len(results)}")
    print(results.head(10).to_string(index=False))


def main():
    parser = argparse.ArgumentParser(
        prog="run_cellpathway",
        description="CellPathway: cell type-specific enhancer enrichment analysis",
    )
    subparsers = parser.add_subparsers(dest="command")

    # --- enrich ---
    enrich = subparsers.add_parser("enrich", help="Run enrichment (Steps 1-7)")
    enrich.add_argument("--enhancer-dir", required=True,
                        help="Directory with <CellType>.hg38.bed files")
    enrich.add_argument("--dnm-file", required=True,
                        help="DNM file (tab-delimited: Chr, Start, End, CADD_PHRED)")
    enrich.add_argument("--output-dir", required=True,
                        help="Output directory")
    enrich.add_argument("--cadd-threshold", type=float, default=10,
                        help="CADD PHRED threshold (default: 10)")
    enrich.add_argument("--bedtools-path", default="bedtools",
                        help="Path to bedtools binary (default: bedtools)")

    # --- tad ---
    tad = subparsers.add_parser("tad", help="TAD annotation (Step 8)")
    tad.add_argument("--overlap-bed", required=True,
                     help="Overlap BED file from enrichment step")
    tad.add_argument("--tad-file", required=True,
                     help="TAD BED file (col4 = 0/1 label)")
    tad.add_argument("--elements-bb", required=True,
                     help="Gene annotation .bb file")
    tad.add_argument("--gene-list", default=None,
                     help="CSV with known disease genes (column 'gene-symbol')")
    tad.add_argument("--output", required=True,
                     help="Output CSV path")

    args = parser.parse_args()
    if args.command is None:
        parser.print_help()
        sys.exit(0)

    if args.command == "enrich":
        cmd_enrich(args)
    elif args.command == "tad":
        cmd_tad(args)


if __name__ == "__main__":
    main()
