#!/usr/bin/env python
"""CellPathway Step 1-7: Cell type-specific enhancer enrichment analysis.

Usage
-----
    python cellpathway_enrich.py \
        --enhancer-dir data/Atlas \
        --dnm-file example/autism_dnm.txt \
        --output-dir example \
        --cadd-threshold 10
"""

import os
import argparse
import subprocess
import pandas as pd
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests


def get_cell_list(enhancer_dir):
    """Get cell type names from enhancer BED file names."""
    return [
        f.replace(".hg38.bed", "")
        for f in os.listdir(enhancer_dir)
        if f.endswith(".hg38.bed")
    ]


def prepare_dnm_bed(dnm_file, output_dir, cadd_threshold):
    """Step 1: Filter DNMs by CADD threshold and write a BED file."""
    gt = pd.read_csv(dnm_file, sep='\t')
    gt = gt[gt['CADD_PHRED'] >= cadd_threshold]
    dnm_bed = gt[['Chr', 'Start', 'End']]
    bed_path = os.path.join(output_dir, f"dnm_{cadd_threshold}.bed")
    dnm_bed.to_csv(bed_path, sep='\t', index=False, header=False)
    print(f"  DNM count (CADD >= {cadd_threshold}): {len(dnm_bed)}")
    return bed_path


def compute_enhancer_lengths(enhancer_dir, final_dir):
    """Step 2: Compute total base-pair length for each enhancer BED file."""
    results = []
    for fname in os.listdir(enhancer_dir):
        if not fname.endswith("hg38.bed"):
            continue
        path = os.path.join(enhancer_dir, fname)
        total = 0
        with open(path) as f:
            for line in f:
                if line.strip() and not line.startswith("#"):
                    parts = line.strip().split()
                    total += int(parts[2]) - int(parts[1])
        clean_name = fname.replace(".hg38.bed", "")
        results.append({"cell": clean_name, "Enhc_bp": total})

    df = pd.DataFrame(results)
    out_csv = os.path.join(final_dir, "Enhc_len.csv")
    df.to_csv(out_csv, index=False, sep='\t')
    return df


def intersect_dnm_enhancers(enhancer_dir, cell_list, dnm_bed_path,
                             overlap_dir, bedtools_path):
    """Step 3: Run bedtools intersect for each cell type."""
    for cell in cell_list:
        cmd = (
            f"{bedtools_path} intersect "
            f"-a {enhancer_dir}/{cell}.hg38.bed "
            f"-b {dnm_bed_path} "
            f"-c > {overlap_dir}/{cell}_dnm.bed"
        )
        try:
            subprocess.run(cmd, shell=True, check=True,
                           capture_output=True, text=True)
        except subprocess.CalledProcessError as e:
            print(f"  [FAIL] {cell}: {e.stderr}")


def count_dnm_overlaps(overlap_dir, final_dir):
    """Step 4: Sum overlap counts per cell type."""
    results = []
    for fname in os.listdir(overlap_dir):
        if not fname.endswith(".bed"):
            continue
        path = os.path.join(overlap_dir, fname)
        try:
            df = pd.read_csv(path, sep='\t', header=None,
                             comment="#", usecols=[3])
            total = df[3].sum()
        except Exception:
            total = 0
        clean = fname.replace("_dnm.bed", "")
        results.append({"cell": clean, "dnm_overlap": total})

    out_df = pd.DataFrame(results)
    out_csv = os.path.join(final_dir, "cell_dnm_overlap.csv")
    out_df.to_csv(out_csv, index=False, sep='\t')
    return out_df


def merge_all_enhancers(enhancer_dir, cell_list, final_dir, bedtools_path):
    """Step 5: Merge all enhancer BED files and return total bp length."""
    bed_files = " ".join(
        f"{enhancer_dir}/{c}.hg38.bed" for c in cell_list
    )
    merged_path = os.path.join(final_dir, "enhc_all_merged.bed")
    cmd = (
        f"cat {bed_files} | sort -k1,1 -k2,2n | "
        f"{bedtools_path} merge -i - > {merged_path}"
    )
    subprocess.run(cmd, shell=True, check=True,
                   capture_output=True, text=True)

    merged = pd.read_csv(merged_path, sep='\t', header=None)
    total_bp = (merged[2] - merged[1]).sum()
    return int(total_bp)


def count_total_dnm_in_enhancers(dnm_bed_path, final_dir, bedtools_path):
    """Step 6: Count total DNMs in all merged enhancer regions."""
    merged_path = os.path.join(final_dir, "enhc_all_merged.bed")
    out_path = os.path.join(final_dir, "all_dnm_enhc.bed")
    cmd = (
        f"{bedtools_path} intersect "
        f"-a {merged_path} -b {dnm_bed_path} -c > {out_path}"
    )
    subprocess.run(cmd, shell=True, check=True,
                   capture_output=True, text=True)
    df = pd.read_csv(out_path, sep='\t', header=None)
    return int(df.iloc[:, 3].sum())


def fisher_enrichment(enhc_len_df, dnm_count_df, total_enhc_bp, total_dnm,
                      final_dir, cadd_threshold):
    """Step 7: Run Fisher's exact test for each cell type and apply FDR."""
    df = pd.merge(enhc_len_df, dnm_count_df, on='cell', how='inner')
    G = total_enhc_bp

    folds, pvals = [], []
    for _, row in df.iterrows():
        A = row["dnm_overlap"]
        B = total_dnm - A
        C = row["Enhc_bp"]
        D = G - C

        fold = (A / C) / (B / D) if B > 0 and D > 0 else 0.0
        table = [[A, C - A], [B, D - B]]
        _, p = fisher_exact(table, alternative="greater")
        folds.append(fold)
        pvals.append(p)

    df["dnm_overlap_other"] = total_dnm - df["dnm_overlap"]
    df["Fold_enrichment"] = folds
    df["P_value"] = pvals
    _, fdr_p, _, _ = multipletests(df["P_value"], method='fdr_bh')
    df["P_FDR"] = fdr_p

    out_path = os.path.join(final_dir, f"enrichment_FC_{cadd_threshold}.csv")
    df.to_csv(out_path, sep=",", index=False)
    print(f"  Results saved to: {out_path}")
    return df


def main():
    parser = argparse.ArgumentParser(
        description="CellPathway: cell type-specific enhancer enrichment analysis (Steps 1-7)",
    )
    parser.add_argument("--enhancer-dir", required=True,
                        help="Directory with <CellType>.hg38.bed files")
    parser.add_argument("--dnm-file", required=True,
                        help="DNM file (tab-delimited: Chr, Start, End, CADD_PHRED)")
    parser.add_argument("--output-dir", required=True,
                        help="Output directory")
    parser.add_argument("--cadd-threshold", type=float, default=10,
                        help="CADD PHRED threshold (default: 10)")
    parser.add_argument("--bedtools-path", default="bedtools",
                        help="Path to bedtools binary (default: bedtools)")
    args = parser.parse_args()

    cadd = args.cadd_threshold
    overlap_dir = os.path.join(args.output_dir, f"dnm_enhc_overlap_cadd_{int(cadd)}")
    final_dir = os.path.join(args.output_dir, f"dnm_enhc_final_cadd_{int(cadd)}")
    os.makedirs(overlap_dir, exist_ok=True)
    os.makedirs(final_dir, exist_ok=True)

    cell_list = get_cell_list(args.enhancer_dir)

    print("Step 1/7: Filtering DNMs by CADD threshold ...")
    dnm_bed = prepare_dnm_bed(args.dnm_file, args.output_dir, cadd)

    print("Step 2/7: Computing enhancer lengths ...")
    enhc_len_df = compute_enhancer_lengths(args.enhancer_dir, final_dir)

    print("Step 3/7: Intersecting DNMs with enhancers (bedtools) ...")
    intersect_dnm_enhancers(args.enhancer_dir, cell_list, dnm_bed,
                            overlap_dir, args.bedtools_path)

    print("Step 4/7: Counting DNM overlaps per cell type ...")
    dnm_count_df = count_dnm_overlaps(overlap_dir, final_dir)

    print("Step 5/7: Merging all enhancer regions ...")
    total_enhc_bp = merge_all_enhancers(args.enhancer_dir, cell_list,
                                         final_dir, args.bedtools_path)

    print("Step 6/7: Counting total DNMs in merged enhancers ...")
    total_dnm = count_total_dnm_in_enhancers(dnm_bed, final_dir,
                                              args.bedtools_path)

    print("Step 7/7: Running Fisher's exact test ...")
    results = fisher_enrichment(enhc_len_df, dnm_count_df, total_enhc_bp,
                                total_dnm, final_dir, cadd)

    print("\nDone! Top 10 enriched cell types:")
    top = results.sort_values("P_FDR").head(10)
    print(top[["cell", "Fold_enrichment", "P_FDR"]].to_string(index=False))


if __name__ == "__main__":
    main()
