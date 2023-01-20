# Assume that current directory contains the DScore.R file and is where the output (category files) should be saved.
# If not, set as required:
setwd(".")

# Load D-Score library. Note: Requires Seurat to be installed.
source("DScore.R")

# This contains a single function: dscore()

# Mode 1: From 10X directory
# Given a 10X path, load the dir with Seurat's Read10X, determine the d-scores and add to the Seurat metadata.

cr_path<-"/path/to/cellranger/output"

# Example:
SeuratObject<-dscore(
  # Required: Path to the cellranger output
  path=cr_path, 
  # Optional. List of features in the 'Antibody Capture' data. Note: Only specify the hashtags used.
  features=c("SBO01", "SBO02", "SBO03", "SBO04", "SBO05"),
  # Optional. Name to save the category assignments to (for loading into Loupe). Leave blank for no category file.
  cat.filename="Category-Cellhashing-Dscore.csv",
  # Optional. Name of the Assignment Category in Loupe
  cat.name="Cell Assignment",
  # Optional. For aggregated Loupe files: increment datasetnumber to append to (an existing) Category file
  # Note: Must be in the same order the 10X samples were aggregated. For 1 sample, leave blank (Default is 1).
  datasetnumber=1
)

# Mode 2: From an existing Seurat Object
# Given a Seurat Object, read the feature barcoding data, determine the d-scores and add to the Seurat metadata.

SeuratObject<-dscore(
  # Required: A Seurat Object
  seurat.obj=SeuratObject,
  # Optional. Set name of the Seurat object's feature barcoding assay. Default is "HTO".
  seurat.fbassay="HTO"
  # Remaining options - features, cat.{filename,name}, datasetnumber - are as above (and optional).
)

## Output from either mode is:
##  1. Seurat Object containing the hash.IDs
##  2. A csv file with assignment categories for loading into Loupe (if cat.filename was specified).

