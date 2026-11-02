# CellPathway: A cell-type-specific enhancer framework for analysis of rare noncoding variants in whole-genome sequencing studies 

CellPathway is an integrative computational framework that combines tissue- and cell-type-specific enhancer annotations with TAD-constrained candidate-gene nomination to facilitate the interpretation of rare noncoding variants and characterize enhancer-enrichment patterns across developmental disorders.


# Installation
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
