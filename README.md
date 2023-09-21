---
title: "tidybulk - part of tidyTranscriptomics"
output: github_document
---

<!-- badges: start -->
[![Lifecycle:maturing](https://img.shields.io/badge/lifecycle-maturing-blue.svg)](https://www.tidyverse.org/lifecycle/#maturing) 
<!-- badges: end -->

# HPCell

## Analysis automated pipeline

Flowchart

https://app.mural.co/t/covid7029/m/covid7029/1656652076667/c47e104697d76b36b8dee3bd05d8de9d96d99efd?sender=udbe1ea99c9618abb07196978

### Clone the repository

```{bash}

git clone git@github.com:susansjy22/HPCell.git

```

### Install Jascap/HPCell package 

```{r, eval=FALSE}

remote::install_github("git@github.com:susansjy22/HPCell.git")

```

### load jascap package 

```{r}

library(jascap)

```

### load input and reference data


```{r}

# Load input data (can be a list of directories or single directory)
library(Seurat)
library(scRNAseq)
file_path = tempfile(tmpdir = "~/HPCell") |> paste0(".rds")
single_cell_data = HeOrganAtlasData(ensembl=FALSE,location=FALSE)[, 1:400]
single_cell_data |> Seurat::as.Seurat(data = NULL) |> saveRDS(file_path)

# Load reference data 
input_reference_path <- "reference_azimuth.rds"
reference_url<- "https://atlas.fredhutch.org/data/nygc/multimodal/pbmc_multimodal.h5seurat"
download.file(reference_url, input_reference_path)
LoadH5Seurat(input_reference_path) |> saveRDS(input_reference_path)
```

### Execute Targets workflow and load results

```{r}
# Provide tissue types: pbmc, solid, atypical
tissue <- "pbmc"

# Create store directory 
store <- "~/HPCell/pipeline_store"

# Default assay name; to be changed to RNA assay
RNA_assay_name <- "originalexp"

# Running the pipeline
preprocessed_seurat = run_targets_pipeline(
    input_data = file_path, 
    store = store,
    input_reference = input_reference_path,
    tissue = tissue,
    computing_resources = crew_controller_local(workers = 1), 
    filtered = TRUE, 
    RNA_assay_name 
)

# Load results
preprocessed_seurat

```
