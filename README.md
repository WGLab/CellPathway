# CellPathway: Cell type-specific enhancer informs pathway-based approaches for association analysis in whole-genome sequencing studies

CellPathway is a framework that directly tests associations between noncoding variants and cell type–specific pathways defined by enhancer activity. 

<img width="710" height="851" alt="image" src="https://github.com/user-attachments/assets/eaefda92-e849-4018-a182-9fe5a3b930fe" />

## Prerequisite

Please refer to `requirements.txt` for required packages.

## Run CellPathway

Please refer to `CellPathway.ipynb` for how to use CellPathway. All inputs and outputs during this example CellPathway run are in the `example` `data` folder.

## Inference
Before running CellPathway, noncoding de novo mutations must be prepared for all samples. We recommend using ANNOVAR to identify and filter rare noncoding variants and to annotate variant pathogenicity using Combined Annotation Dependent Depletion (CADD) scores. These annotated variants serve as the input for downstream cell type–specific pathway analysis in CellPathway.

### Install ANNOVAR
Typically you will go to the [ANNOVAR website](https://annovar.openbioinformatics.org/en/latest/), fill in a registration form, and download the package there. When you have requested the ANNOVAR from the website and downloaded it, you will receive a compressed file ```annovar.latest.tar.gz```, you will need to unzip it. Then follow the user guide to install ANNOVAR.

### Install CADD
Please follow [CADD](https://github.com/kircherlab/CADD-scripts) repository for instructions on how to install and run CADD.
