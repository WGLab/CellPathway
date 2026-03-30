# CellPathway: Cell type-specific enhancer informs pathway-based approaches for association analysis in whole-genome sequencing studies

CellPathway is a framework that directly tests associations between noncoding variants and cell type-specific pathways defined by enhancer activity.

<img width="710" height="851" alt="image" src="https://github.com/user-attachments/assets/eaefda92-e849-4018-a182-9fe5a3b930fe" />

## Installation

### Step 1: Clone the repository

```bash
git clone https://github.com/WGLab/CellPathway.git
cd CellPathway
```

### Step 2: Create a conda environment

```bash
conda create -n cellpathway python=3.10 -y
conda activate cellpathway
```

### Step 3: Install bedtools

```bash
conda install -c bioconda bedtools -y
```

### Step 4: Install Python dependencies

```bash
pip install -r requirements.txt
```

## Quick start

### Step 1-7: Run enrichment analysis

```bash
python run_cellpathway.py enrich \
    --enhancer-dir data/Atlas \
    --dnm-file example/autism_dnm.txt \
    --output-dir example \
    --cadd-threshold 10
```

### Step 8: TAD annotation

After enrichment, run TAD annotation for a specific cell type using the overlap BED file generated in the previous step:

```bash
python run_cellpathway.py tad \
    --overlap-bed example/dnm_enhc_overlap_cadd_10/Fetal_brain_dnm.bed \
    --tad-file data/tad_w_boundary_08.bed \
    --elements-bb data/genes_w_noncoding.bb \
    --gene-list example/SFARI_Gene.csv \
    --output example/tad_Fetal_brain_autism.csv
```

The `--overlap-bed` path comes from the enrichment output. Replace `Fetal_brain` with any cell type from your results.

## Input data

| File | Description |
|------|-------------|
| `example/autism_dnm.txt` | De novo mutations with CADD scores (tab-delimited: Chr, Start, End, CADD_PHRED) |
| `data/Atlas/` | Cell type-specific enhancer BED files (51 tissues, `*.hg38.bed`) |
| `data/Single_cell/` | Single-cell enhancer BED files (169 cell types, `*.hg38.bed`) |
| `data/tad_w_boundary_08.bed` | TAD regions and boundaries (hg38) |
| `data/genes_w_noncoding.bb` | Gene annotations in BigWig format |
| `example/SFARI_Gene.csv` | SFARI autism-associated gene list |

## Output

- **Enrichment table** (`enrichment_FC_<cadd>.csv`): cell type, enhancer bp, DNM overlap count, fold enrichment, p-value, FDR-adjusted p-value
- **TAD annotation** (`tad_<cell>_autism.csv`): enhancer-to-TAD mapping with gene names and known disease gene overlaps

## Preprocessing

Before running CellPathway, noncoding de novo mutations must be prepared for all samples. We recommend using ANNOVAR to identify and filter rare noncoding variants and to annotate variant pathogenicity using Combined Annotation Dependent Depletion (CADD) scores.

- [ANNOVAR](https://annovar.openbioinformatics.org/en/latest/)
- [CADD](https://github.com/kircherlab/CADD-scripts)

## License

MIT
