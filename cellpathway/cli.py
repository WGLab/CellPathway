"""Command-line interface for CellPathway."""

import argparse
import sys

from cellpathway.core import CellPathway
from cellpathway.tad import annotate_tad


def main():
    parser = argparse.ArgumentParser(
        prog="cellpathway",
        description="CellPathway: cell type-specific enhancer enrichment analysis",
    )
    subparsers = parser.add_subparsers(dest="command")

    # --- enrich sub-command ---
    enrich = subparsers.add_parser(
        "enrich",
        help="Run the enrichment pipeline (Steps 1-7)",
    )
    enrich.add_argument("--enhancer-dir", required=True,
                        help="Directory with <CellType>.hg38.bed files")
    enrich.add_argument("--dnm-file", required=True,
                        help="Tab-delimited DNM file (Chr, Start, End, CADD_PHRED)")
    enrich.add_argument("--output-dir", required=True,
                        help="Output directory")
    enrich.add_argument("--cadd-threshold", type=float, default=10,
                        help="CADD PHRED score threshold (default: 10)")
    enrich.add_argument("--bedtools-path", default="bedtools",
                        help="Path to bedtools binary (default: bedtools)")

    # --- tad sub-command ---
    tad = subparsers.add_parser(
        "tad",
        help="Annotate enriched enhancers with TAD regions and genes (Step 8)",
    )
    tad.add_argument("--overlap-bed", required=True,
                     help="Bedtools overlap BED file for one cell type")
    tad.add_argument("--tad-file", required=True,
                     help="TAD BED file (4 columns; col4 = 0/1 label)")
    tad.add_argument("--elements-bb", required=True,
                     help="Gene annotation .bb file")
    tad.add_argument("--gene-list", default=None,
                     help="CSV file with known disease genes "
                          "(column 'gene-symbol')")
    tad.add_argument("--output", required=True,
                     help="Output CSV path")

    args = parser.parse_args()

    if args.command is None:
        parser.print_help()
        sys.exit(0)

    if args.command == "enrich":
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

    elif args.command == "tad":
        import pandas as pd
        gene_list = None
        if args.gene_list:
            genes_df = pd.read_csv(args.gene_list)
            gene_list = genes_df["gene-symbol"].tolist()

        annotate_tad(
            overlap_bed_path=args.overlap_bed,
            tad_path=args.tad_file,
            elements_bb_path=args.elements_bb,
            gene_list=gene_list,
            output_path=args.output,
        )


if __name__ == "__main__":
    main()
