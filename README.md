tidybulk - part of tidyTranscriptomics
================

<!-- badges: start -->

[![Lifecycle:maturing](https://img.shields.io/badge/lifecycle-maturing-blue.svg)](https://www.tidyverse.org/lifecycle/#maturing)
<!-- badges: end -->

# HPCell

## Analysis automated pipeline

Flowchart

<https://app.mural.co/t/covid7029/m/covid7029/1656652076667/c47e104697d76b36b8dee3bd05d8de9d96d99efd?sender=udbe1ea99c9618abb07196978>

Install Jascap/HPCell package

``` r
remote::install_github("stemangiola/HPCell")
```

load jascap package

``` r
library(HPCell)
```

load input and reference data

``` r
# Load input data (can be a list of directories or single directory)
library(Seurat)
library(scRNAseq)
input_data_path =  tempfile(tmpdir = ".") |> paste0(".rds")
HeOrganAtlasData(ensembl=FALSE,location=FALSE)[, 1:400] |>
  as.Seurat(data = NULL) |>
  saveRDS(input_data_path)
```

Execute Targets workflow and load results

``` r
# Running the pipeline
preprocessed_seurat = run_targets_pipeline(
    input_data = input_data_path,
    tissue = "pbmc",
    filter_input = TRUE,
    RNA_assay_name = "originalexp",
    sample_column = "Tissue", 
    debug_step = "annotation_label_transfer_tbl"
)
```

    ## Warning: Targets and globals must have unique names. Ignoring global objects
    ## that conflict with target names: sample_column, tissue, filter_input. Warnings
    ## like this one are important, but if you must suppress them, you can do so with
    ## Sys.setenv(TAR_WARN = "false").

    ## ✔ skip target reference_file

    ## ✔ skip target tissue_file
    ## ✔ skip target sample_column_file
    ## ✔ skip target sample_column
    ## ✔ skip target tissue
    ## ✔ skip target reference_label_fine
    ## ✔ skip target read_file
    ## ✔ skip branch input_read_46201ef3
    ## ✔ skip pattern input_read
    ## ✔ skip branch input_read_RNA_assay_363b64d4
    ## ✔ skip pattern input_read_RNA_assay
    ## ✔ skip target file
    ## ✔ skip target reference_read
    ## ✔ skip target reference_label_coarse
    ## ✔ skip target filtered_file
    ## ✔ skip target filter_input
    ## ✔ skip branch empty_droplets_tbl_58837869
    ## ✔ skip pattern empty_droplets_tbl
    ## ✔ skip branch annotation_label_transfer_tbl_f47e871e
    ## ✔ skip pattern annotation_label_transfer_tbl
    ## ✔ skip branch alive_identification_tbl_31879d89
    ## ✔ skip pattern alive_identification_tbl
    ## ✔ skip branch doublet_identification_tbl_4885fe19
    ## ✔ skip pattern doublet_identification_tbl
    ## ✔ skip branch cell_cycle_score_tbl_f47e871e
    ## ✔ skip pattern cell_cycle_score_tbl
    ## ✔ skip branch non_batch_variation_removal_S_6c860224
    ## ✔ skip pattern non_batch_variation_removal_S
    ## ✔ skip branch preprocessing_output_S_2745fad8
    ## ✔ skip pattern preprocessing_output_S
    ## ✔ skip target pseudobulk_preprocessing_SE

    ## ✔ skip pipeline [0.113 seconds]

    ## HPCell says: you can read your output executing tar_read(preprocessing_output_S, store = "./")

``` r
# Load results
preprocessed_seurat
```

    ## $preprocessing_output_S_2745fad8
    ## # A Seurat-tibble abstraction: 284 × 46
    ## # [90mFeatures=9552 | Cells=284 | Active assay=SCT | Assays=RNA, SCT[0m
    ##    .cell    orig.ident nCount_originalexp nFeature_originalexp Tissue nCount_RNA
    ##    <chr>    <fct>                   <dbl>                <int> <chr>       <dbl>
    ##  1 Bladder… Bladder                  1133                  592 Bladd…       1133
    ##  2 Bladder… Bladder                  3495                 1374 Bladd…       3495
    ##  3 Bladder… Bladder                  1297                  797 Bladd…       1297
    ##  4 Bladder… Bladder                  2071                  953 Bladd…       2071
    ##  5 Bladder… Bladder                  2166                 1102 Bladd…       2166
    ##  6 Bladder… Bladder                  2486                 1185 Bladd…       2486
    ##  7 Bladder… Bladder                  3175                 1418 Bladd…       3175
    ##  8 Bladder… Bladder                  1847                  932 Bladd…       1847
    ##  9 Bladder… Bladder                  2546                 1104 Bladd…       2546
    ## 10 Bladder… Bladder                   969                  574 Bladd…        969
    ## # ℹ 274 more rows
    ## # ℹ 40 more variables: nFeature_RNA <int>, percent.mito <dbl>,
    ## #   RNA_snn_res.orig <int>, seurat_clusters <int>,
    ## #   Cell_type_in_each_tissue <chr>, Cell_type_in_merged_data <chr>,
    ## #   reclustered.broad <chr>, reclustered.fine <chr>, Total <int>,
    ## #   LogProb <dbl>, PValue <dbl>, Limited <lgl>, FDR <dbl>, empty_droplet <lgl>,
    ## #   rank <dbl>, total <int>, fitted <dbl>, knee <dbl>, inflection <dbl>, …

Include reference dataset for azimuth annotation

``` r
# Load reference data
input_reference_path <- "reference_azimuth.rds"
reference_url<- "https://atlas.fredhutch.org/data/nygc/multimodal/pbmc_multimodal.h5seurat"
download.file(reference_url, input_reference_path)
LoadH5Seurat(input_reference_path) |> saveRDS(input_reference_path)
```
