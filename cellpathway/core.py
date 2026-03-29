"""Core CellPathway pipeline."""

import os
import subprocess
import pandas as pd
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests
from pybedtools import BedTool


class CellPathway:
    """Cell type-specific enhancer enrichment analysis pipeline.

    Parameters
    ----------
    enhancer_dir : str
        Directory containing cell type-specific enhancer BED files
        (named ``<CellType>.hg38.bed``).
    dnm_file : str
        Path to de novo mutation file (tab-delimited with columns:
        Chr, Start, End, and CADD_PHRED).
    output_dir : str
        Directory for output files.
    cadd_threshold : int or float
        CADD PHRED score threshold for filtering variants (default 10).
    bedtools_path : str
        Path to bedtools binary. Defaults to ``"bedtools"`` (i.e. on PATH).
    """

    def __init__(self, enhancer_dir, dnm_file, output_dir, cadd_threshold=10,
                 bedtools_path="bedtools"):
        self.enhancer_dir = enhancer_dir
        self.dnm_file = dnm_file
        self.output_dir = output_dir
        self.cadd_threshold = cadd_threshold
        self.bedtools_path = bedtools_path

        self.overlap_dir = os.path.join(output_dir, f"dnm_enhc_overlap_cadd_{cadd_threshold}")
        self.final_dir = os.path.join(output_dir, f"dnm_enhc_final_cadd_{cadd_threshold}")
        os.makedirs(self.overlap_dir, exist_ok=True)
        os.makedirs(self.final_dir, exist_ok=True)

        self._cell_list = None

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def run(self):
        """Run the full enrichment pipeline and return the results DataFrame."""
        print("Step 1/7: Filtering DNMs by CADD threshold ...")
        dnm_bed = self.prepare_dnm_bed()

        print("Step 2/7: Computing enhancer lengths ...")
        enhc_len_df = self.compute_enhancer_lengths()

        print("Step 3/7: Intersecting DNMs with enhancers (bedtools) ...")
        self.intersect_dnm_enhancers(dnm_bed)

        print("Step 4/7: Counting DNM overlaps per cell type ...")
        dnm_count_df = self.count_dnm_overlaps()

        print("Step 5/7: Merging all enhancer regions ...")
        total_enhc_bp = self.merge_all_enhancers()

        print("Step 6/7: Counting total DNMs in merged enhancers ...")
        total_dnm = self.count_total_dnm_in_enhancers(dnm_bed)

        print("Step 7/7: Running Fisher's exact test ...")
        results = self.fisher_enrichment(enhc_len_df, dnm_count_df,
                                         total_enhc_bp, total_dnm)
        print("Done.")
        return results

    def prepare_dnm_bed(self):
        """Filter DNMs by CADD threshold and write a BED file.

        Returns the path to the generated BED file.
        """
        gt = pd.read_csv(self.dnm_file, sep='\t')
        gt = gt[gt['CADD_PHRED'] >= self.cadd_threshold]
        dnm_bed = gt[['Chr', 'Start', 'End']]
        bed_path = os.path.join(self.output_dir,
                                f"dnm_{self.cadd_threshold}.bed")
        dnm_bed.to_csv(bed_path, sep='\t', index=False, header=False)
        print(f"  DNM count (CADD >= {self.cadd_threshold}): {len(dnm_bed)}")
        return bed_path

    def compute_enhancer_lengths(self):
        """Compute total base-pair length for each enhancer BED file."""
        results = []
        for fname in os.listdir(self.enhancer_dir):
            if not fname.endswith("hg38.bed"):
                continue
            path = os.path.join(self.enhancer_dir, fname)
            total = 0
            with open(path) as f:
                for line in f:
                    if line.strip() and not line.startswith("#"):
                        parts = line.strip().split()
                        total += int(parts[2]) - int(parts[1])
            clean_name = fname.replace(".hg38.bed", "")
            results.append({"cell": clean_name, "Enhc_bp": total})

        df = pd.DataFrame(results)
        out_csv = os.path.join(self.final_dir, "Enhc_len.csv")
        df.to_csv(out_csv, index=False, sep='\t')
        return df

    @property
    def cell_list(self):
        """List of cell type names derived from enhancer BED file names."""
        if self._cell_list is None:
            self._cell_list = [
                f.replace(".hg38.bed", "")
                for f in os.listdir(self.enhancer_dir)
                if f.endswith(".hg38.bed")
            ]
        return self._cell_list

    def intersect_dnm_enhancers(self, dnm_bed_path):
        """Run ``bedtools intersect`` for each cell type."""
        for cell in self.cell_list:
            cmd = (
                f"{self.bedtools_path} intersect "
                f"-a {self.enhancer_dir}/{cell}.hg38.bed "
                f"-b {dnm_bed_path} "
                f"-c > {self.overlap_dir}/{cell}_dnm.bed"
            )
            try:
                subprocess.run(cmd, shell=True, check=True,
                               capture_output=True, text=True)
            except subprocess.CalledProcessError as e:
                print(f"  [FAIL] {cell}: {e.stderr}")

    def count_dnm_overlaps(self):
        """Sum overlap counts per cell type from bedtools output."""
        results = []
        for fname in os.listdir(self.overlap_dir):
            if not fname.endswith(".bed"):
                continue
            path = os.path.join(self.overlap_dir, fname)
            try:
                df = pd.read_csv(path, sep='\t', header=None,
                                 comment="#", usecols=[3])
                total = df[3].sum()
            except Exception:
                total = 0
            clean = fname.replace("_dnm.bed", "")
            results.append({"cell": clean, "dnm_overlap": total})

        out_df = pd.DataFrame(results)
        out_csv = os.path.join(self.final_dir, "cell_dnm_overlap.csv")
        out_df.to_csv(out_csv, index=False, sep='\t')
        return out_df

    def merge_all_enhancers(self):
        """Merge all enhancer BED files and return total base-pair length."""
        bed_files = " ".join(
            f"{self.enhancer_dir}/{c}.hg38.bed" for c in self.cell_list
        )
        merged_path = os.path.join(self.final_dir, "enhc_all_merged.bed")
        cmd = (
            f"cat {bed_files} | sort -k1,1 -k2,2n | "
            f"{self.bedtools_path} merge -i - > {merged_path}"
        )
        subprocess.run(cmd, shell=True, check=True,
                       capture_output=True, text=True)

        merged = pd.read_csv(merged_path, sep='\t', header=None)
        total_bp = (merged[2] - merged[1]).sum()
        return int(total_bp)

    def count_total_dnm_in_enhancers(self, dnm_bed_path):
        """Count total DNMs falling in any merged enhancer region."""
        merged_path = os.path.join(self.final_dir, "enhc_all_merged.bed")
        out_path = os.path.join(self.final_dir, "all_dnm_enhc.bed")
        cmd = (
            f"{self.bedtools_path} intersect "
            f"-a {merged_path} -b {dnm_bed_path} -c > {out_path}"
        )
        subprocess.run(cmd, shell=True, check=True,
                       capture_output=True, text=True)
        df = pd.read_csv(out_path, sep='\t', header=None)
        return int(df.iloc[:, 3].sum())

    def fisher_enrichment(self, enhc_len_df, dnm_count_df,
                          total_enhc_bp, total_dnm):
        """Run Fisher's exact test for each cell type and apply FDR."""
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

        out_path = os.path.join(
            self.final_dir, f"enrichment_FC_{self.cadd_threshold}.csv"
        )
        df.to_csv(out_path, sep=",", index=False)
        print(f"  Results saved to: {out_path}")
        return df
