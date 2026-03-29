# CellPathway: Cell type-specific enhancer informs pathway-based approaches for association analysis in whole-genome sequencing studies

CellPathway is a framework that directly tests associations between noncoding variants and cell type-specific pathways defined by enhancer activity.

<img width="710" height="851" alt="image" src="https://github.com/user-attachments/assets/eaefda92-e849-4018-a182-9fe5a3b930fe" />

## Installation

### Prerequisites

- Python >= 3.9
- [bedtools](https://bedtools.readthedocs.io/) must be installed and available on your `PATH`

### Install from source

```bash
git clone https://github.com/<your-username>/CellPathway.git
cd CellPathway
pip install .
```

For development (editable install):

```bash
pip install -e .
```

## Quick start

### Python API

```python
from cellpathway import CellPathway
from cellpathway.tad import annotate_tad

# Step 1-7: Enrichment analysis
cp = CellPathway(
    enhancer_dir="data/Atlas",
    dnm_file="example/autism_dnm.txt",
    output_dir="example",
    cadd_threshold=10,
)
results = cp.run()

# Show top enriched cell types
print(results.sort_values("P_FDR").head(10))

# Step 8: TAD annotation for a specific cell type
import pandas as pd
sfari = pd.read_csv("example/SFARI_Gene.csv")

tad_results = annotate_tad(
    overlap_bed_path="example/dnm_enhc_overlap_cadd_10/Fetal_brain_dnm.bed",
    tad_path="data/tad_w_boundary_08.bed",
    elements_bb_path="data/genes_w_noncoding.bb",
    gene_list=sfari["gene-symbol"].tolist(),
    output_path="example/tad_brain_autism.csv",
)
```

### Command line

```bash
# Run enrichment (Steps 1-7)
cellpathway enrich \
    --enhancer-dir data/Atlas \
    --dnm-file example/autism_dnm.txt \
    --output-dir example \
    --cadd-threshold 10

# TAD annotation (Step 8)
cellpathway tad \
    --overlap-bed example/dnm_enhc_overlap_cadd_10/Fetal_brain_dnm.bed \
    --tad-file data/tad_w_boundary_08.bed \
    --elements-bb data/genes_w_noncoding.bb \
    --gene-list example/SFARI_Gene.csv \
    --output example/tad_brain_autism.csv
```

### Example script

A complete example using the bundled autism dataset:

```bash
python run_cellpathway.py
```

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
