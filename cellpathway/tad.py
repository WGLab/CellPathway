"""TAD annotation and gene mapping."""

import os
import pandas as pd
from pybedtools import BedTool
import pyBigWig


def superset_tad(query_region, tad_r, tad_b):
    """Determine whether a query region is within a TAD or crosses a boundary.

    Parameters
    ----------
    query_region : list
        ``[chrom, start, end]``
    tad_r : pd.DataFrame
        TAD regions (label == 0).
    tad_b : pd.DataFrame
        TAD boundaries (label == 1).

    Returns
    -------
    region : list
        ``[chrom, start, end]`` of the covering TAD.
    status : str
        One of ``"in_TAD"``, ``"cross_TAD_boundary"``, or ``"outside_TAD"``.
    """
    if isinstance(query_region, (list, tuple)):
        query_region = ' '.join(str(i) for i in query_region)
    query_bed = BedTool(query_region, from_string=True)
    tad_r_bed = BedTool.from_dataframe(tad_r)
    tad_b_bed = BedTool.from_dataframe(tad_b)

    inter_b = tad_b_bed.intersect(query_bed, wa=True, wb=True).to_dataframe()

    if inter_b.shape[0] == 0:
        inter_r = tad_r_bed.intersect(query_bed, wa=True).to_dataframe()
        if inter_r.shape[0] == 0:
            return [None, None, None], "outside_TAD"
        return inter_r.values.tolist()[0], "in_TAD"

    boundary_selected = ' '.join([
        inter_b.chrom.iloc[0],
        str(min(inter_b.start)),
        str(max(inter_b.end))
    ])
    boundary_bed = BedTool(boundary_selected, from_string=True)
    fu = boundary_bed.closest(tad_r_bed, D='ref', fu=True).to_dataframe()
    fd = boundary_bed.closest(tad_r_bed, D='ref', fd=True).to_dataframe()
    region = [fu.values[0][0], fu.values[0][4], fd.values[0][5]]
    return region, "cross_TAD_boundary"


def tad_elements(element_bb, tad_region):
    """Get all genes/elements within a TAD region.

    Parameters
    ----------
    element_bb : pyBigWig file handle
        Opened ``.bb`` file with gene annotations.
    tad_region : list
        ``[chrom, start, end]``

    Returns
    -------
    gene_names : list of str
        Gene names found in the TAD.
    elements_df : pd.DataFrame or None
        All elements in the region.
    """
    if None in tad_region:
        return [], None

    try:
        entries = element_bb.entries(
            tad_region[0], int(tad_region[1]), int(tad_region[2])
        )
    except Exception:
        return [], None

    gene_names = []
    rows = []
    for entry in entries:
        name, indicator = entry[2].split('\t')
        rows.append([tad_region[0], entry[0], entry[1], name])
        if indicator == '1':
            gene_names.append(name)

    elements_df = pd.DataFrame(rows, columns=["chrom", "start", "end", "feature"])
    return gene_names, elements_df


def annotate_tad(overlap_bed_path, tad_path, elements_bb_path,
                 gene_list=None, output_path=None):
    """Map enhancer-DNM overlaps to TADs and identify genes.

    Parameters
    ----------
    overlap_bed_path : str
        Path to a bedtools overlap BED file (4-column; col 4 = overlap count).
    tad_path : str
        Path to TAD BED file (4-column; col 4 = label 0/1).
    elements_bb_path : str
        Path to ``.bb`` gene annotation file.
    gene_list : list of str, optional
        Known disease gene list for overlap (e.g. SFARI genes).
    output_path : str, optional
        If provided, save results to this CSV path.

    Returns
    -------
    pd.DataFrame
        Annotated results with TAD regions and gene names.
    """
    # Load TAD data
    tad = pd.read_csv(tad_path, sep="\t", header=None,
                       names=["chrom", "start", "end", "label"])
    tad_r = tad.loc[tad["label"].astype(str) == "0",
                     ["chrom", "start", "end"]].reset_index(drop=True)
    tad_b = tad.loc[tad["label"].astype(str) == "1",
                     ["chrom", "start", "end"]].reset_index(drop=True)

    element_bb = pyBigWig.open(elements_bb_path, "r")

    # Load overlap BED and filter to regions with overlaps
    overlap = pd.read_csv(overlap_bed_path, sep="\t", header=None)
    enhancers = overlap[overlap[3] > 0][[0, 1, 2]]
    enhancers.columns = ["chrom", "start", "end"]

    gene_set = set(g.upper() for g in gene_list) if gene_list else None

    results = []
    for _, row in enhancers.iterrows():
        enhc_tri = [row["chrom"], int(row["start"]), int(row["end"])]
        tad_region, status = superset_tad(enhc_tri, tad_r, tad_b)
        gene_names, _ = tad_elements(element_bb, tad_region)

        record = {
            "enh_chr": row["chrom"],
            "enh_start": row["start"],
            "enh_end": row["end"],
            "tad_chr": tad_region[0],
            "tad_start": tad_region[1],
            "tad_end": tad_region[2],
            "status": status,
            "gene_names": ",".join(gene_names) if gene_names else None,
            "n_genes": len(gene_names),
        }

        if gene_set is not None:
            genes_upper = [g.strip().upper() for g in gene_names]
            overlap_genes = [g for g in genes_upper if g in gene_set]
            record["known_gene_overlap"] = ",".join(overlap_genes)

        results.append(record)

    element_bb.close()
    results_df = pd.DataFrame(results)

    if output_path:
        results_df.to_csv(output_path, index=False)
        print(f"TAD annotation saved to: {output_path}")

    return results_df
