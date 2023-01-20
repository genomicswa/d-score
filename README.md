# dscore :: a deconvolution method for single-cell RNA-Seq

This repo contains the R code to deconvolute single cells tagged with SBOs back to their original assignment.

## Usage

D-Score requires Seurat to run, and has been developed with Seurat v4.0.2.

 - `DScore.R`: Contains the dscore function
 - `Example.R`: Sample R code on how to use D-Score for SBO deconvolution

### Mode 1: from cellranger output

A simple example:
```R
source(DScore.R)

SeuratObject<-dscore(
  path="/path/to/cellranger/output", 
  # List of features to include
  features=c("SBO01", "SBO02", "SBO03", "SBO04", "SBO05"),
  # Name of the assignment category csv output file
  cat.filename="Category-Cellhashing-Dscore.csv"
)
```

### Mode 2: from existing Seurat object

D-Score can also be run from an existing Seurat object (see `Example.R` for details). To create Seurat Objects and load the HTO assay, you may find [this dehashing tutorial](https://satijalab.org/seurat/articles/hashing_vignette.html) useful, particularly the section on "Adding HTO data as an independent assay".

## Output

 1. `SeuratObject` with 'hash.ID' in the metadata
 2. (if specified) CSV file with hash assignment information for loading into Loupe Cell Browser

## References

*Publication in progress..*
