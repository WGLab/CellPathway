# CellPathway: A cell-type-specific enhancer framework for analysis of rare noncoding variants in whole-genome sequencing studies 

CellPathway is an integrative computational framework that combines tissue- and cell-type-specific enhancer annotations with TAD-constrained candidate-gene nomination to facilitate the interpretation of rare noncoding variants and characterize enhancer-enrichment patterns across developmental disorders.


## Installation
We recommend using Conda to set up the environment. If Conda is not installed, run the following commands in Linux to install it.

```
curl -O https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
```

After conda is installed successfully, CellPathway sources can be downloaded:
```
git clone https://github.com/WGLab/CellPathway.git
cd cellpathway
conda create -n cellpathway python=3.10
conda activate cellpathway
conda install -c conda-forge -c bioconda bedtools numpy pandas scipy statsmodels pybigwig -y
```
## Quick start

### Run enrichment analysis

### Input data
- Enhancer and TAD coordinates are GRCh38 BED coordinates.
- Custom DNM inputs must be tab-delimited BED files with at least three columns and chromosome names matching the enhancer files (chr1, chr2, ...).
- To analyze a custom DNM file, pass --dnv file.

Check resource
```
python scripts/run_enrichment.py --resource atlas
python scripts/run_enrichment.py --resource scenhancer
```
Run example
```
python scripts/run_enrichment.py \
  --resource atlas \
  --dnv example/example.bed
```
## TAD-constrained candidate-gene annotation
After selecting a phenotype-relevant enriched enhancer set, annotate its overlapping DNMs. For example, fetal-brain enhancers:
```
python scripts/annotate_tad.py \
  --dnv example/example.bed \
  --resource atlas \
  --enhancer-set Fetal_brain \
  --tad brain \
  --output results/tad/example_fetal_brain.tsv
```
For craniofacial analyses, use --tad cncc. The output contains the overlapping enhancer interval, TAD interval, boundary-crossing status, and all protein-coding genes in that TAD.

TAD co-membership nominates genes in the same regulatory neighborhood; it does not establish a direct enhancer-target relationship.

## Preprocessing

Before running CellPathway, noncoding de novo mutations must be prepared for all samples. We recommend using ANNOVAR to identify and filter rare noncoding variants and to annotate variant pathogenicity using Combined Annotation Dependent Depletion (CADD) scores.

- [ANNOVAR](https://annovar.openbioinformatics.org/en/latest/)
- [CADD](https://github.com/kircherlab/CADD-scripts)

