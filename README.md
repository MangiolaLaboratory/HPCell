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

    ## The legacy packages maptools, rgdal, and rgeos, underpinning the sp package,
    ## which was just loaded, will retire in October 2023.
    ## Please refer to R-spatial evolution reports for details, especially
    ## https://r-spatial.org/r/2023/05/15/evolution4.html.
    ## It may be desirable to make the sf package available;
    ## package maintainers should consider adding sf to Suggests:.
    ## The sp package is now running under evolution status 2
    ##      (status 2 uses the sf package in place of rgdal)

    ## Warning: replacing previous import 'tidySingleCellExperiment::plot_ly' by
    ## 'tidySummarizedExperiment::plot_ly' when loading 'HPCell'

    ## Warning: replacing previous import 'tidySingleCellExperiment::tidy' by
    ## 'tidySummarizedExperiment::tidy' when loading 'HPCell'

    ## Warning: replacing previous import 'tidySingleCellExperiment::bind_cols' by
    ## 'tidySummarizedExperiment::bind_cols' when loading 'HPCell'

    ## Warning: replacing previous import 'tidySingleCellExperiment::bind_rows' by
    ## 'tidySummarizedExperiment::bind_rows' when loading 'HPCell'

    ## Warning: replacing previous import 'tidySingleCellExperiment::count' by
    ## 'tidySummarizedExperiment::count' when loading 'HPCell'

    ## Warning: replacing previous import 'tidySummarizedExperiment::tidy' by
    ## 'tidyseurat::tidy' when loading 'HPCell'

    ## Warning: replacing previous import 'tidySummarizedExperiment::plot_ly' by
    ## 'tidyseurat::plot_ly' when loading 'HPCell'

    ## Warning: replacing previous import 'tidySingleCellExperiment::aggregate_cells'
    ## by 'tidyseurat::aggregate_cells' when loading 'HPCell'

    ## Warning: replacing previous import 'tidySingleCellExperiment::join_transcripts'
    ## by 'tidyseurat::join_transcripts' when loading 'HPCell'

    ## Warning: replacing previous import 'tidySummarizedExperiment::count' by
    ## 'dplyr::count' when loading 'HPCell'

load input and reference data

``` r
# Load input data (can be a list of directories or single directory)
library(Seurat)
```

    ## Attaching SeuratObject

    ## 'SeuratObject' was built under R 4.3.0 but the current version is
    ## 4.3.1; it is recomended that you reinstall 'SeuratObject' as the ABI
    ## for R may have changed

``` r
library(scRNAseq)
```

    ## Loading required package: SingleCellExperiment

    ## Loading required package: SummarizedExperiment

    ## Loading required package: MatrixGenerics

    ## Loading required package: matrixStats

    ## 
    ## Attaching package: 'matrixStats'

    ## The following objects are masked from 'package:Biobase':
    ## 
    ##     anyMissing, rowMedians

    ## 
    ## Attaching package: 'MatrixGenerics'

    ## The following objects are masked from 'package:matrixStats':
    ## 
    ##     colAlls, colAnyNAs, colAnys, colAvgsPerRowSet, colCollapse,
    ##     colCounts, colCummaxs, colCummins, colCumprods, colCumsums,
    ##     colDiffs, colIQRDiffs, colIQRs, colLogSumExps, colMadDiffs,
    ##     colMads, colMaxs, colMeans2, colMedians, colMins, colOrderStats,
    ##     colProds, colQuantiles, colRanges, colRanks, colSdDiffs, colSds,
    ##     colSums2, colTabulates, colVarDiffs, colVars, colWeightedMads,
    ##     colWeightedMeans, colWeightedMedians, colWeightedSds,
    ##     colWeightedVars, rowAlls, rowAnyNAs, rowAnys, rowAvgsPerColSet,
    ##     rowCollapse, rowCounts, rowCummaxs, rowCummins, rowCumprods,
    ##     rowCumsums, rowDiffs, rowIQRDiffs, rowIQRs, rowLogSumExps,
    ##     rowMadDiffs, rowMads, rowMaxs, rowMeans2, rowMedians, rowMins,
    ##     rowOrderStats, rowProds, rowQuantiles, rowRanges, rowRanks,
    ##     rowSdDiffs, rowSds, rowSums2, rowTabulates, rowVarDiffs, rowVars,
    ##     rowWeightedMads, rowWeightedMeans, rowWeightedMedians,
    ##     rowWeightedSds, rowWeightedVars

    ## The following object is masked from 'package:Biobase':
    ## 
    ##     rowMedians

    ## Loading required package: GenomicRanges

    ## Loading required package: stats4

    ## Loading required package: S4Vectors

    ## 
    ## Attaching package: 'S4Vectors'

    ## The following object is masked from 'package:utils':
    ## 
    ##     findMatches

    ## The following objects are masked from 'package:base':
    ## 
    ##     expand.grid, I, unname

    ## Loading required package: IRanges

    ## Loading required package: GenomeInfoDb

    ## 
    ## Attaching package: 'SummarizedExperiment'

    ## The following object is masked from 'package:SeuratObject':
    ## 
    ##     Assays

    ## The following object is masked from 'package:Seurat':
    ## 
    ##     Assays

``` r
input_data_path =  tempfile(tmpdir = "~") |> paste0(".rds")
HeOrganAtlasData(ensembl=FALSE,location=FALSE)[, 1:400] |> 
  as.Seurat(data = NULL) |> 
  saveRDS(input_data_path) 
```

    ## see ?scRNAseq and browseVignettes('scRNAseq') for documentation

    ## loading from cache

    ## see ?scRNAseq and browseVignettes('scRNAseq') for documentation

    ## loading from cache

    ## see ?scRNAseq and browseVignettes('scRNAseq') for documentation

    ## loading from cache

    ## see ?scRNAseq and browseVignettes('scRNAseq') for documentation

    ## loading from cache

    ## see ?scRNAseq and browseVignettes('scRNAseq') for documentation

    ## loading from cache

    ## see ?scRNAseq and browseVignettes('scRNAseq') for documentation

    ## loading from cache

    ## see ?scRNAseq and browseVignettes('scRNAseq') for documentation

    ## loading from cache

    ## see ?scRNAseq and browseVignettes('scRNAseq') for documentation

    ## loading from cache

    ## see ?scRNAseq and browseVignettes('scRNAseq') for documentation

    ## loading from cache

    ## see ?scRNAseq and browseVignettes('scRNAseq') for documentation

    ## loading from cache

    ## see ?scRNAseq and browseVignettes('scRNAseq') for documentation

    ## loading from cache

    ## see ?scRNAseq and browseVignettes('scRNAseq') for documentation

    ## loading from cache

    ## see ?scRNAseq and browseVignettes('scRNAseq') for documentation

    ## loading from cache

    ## see ?scRNAseq and browseVignettes('scRNAseq') for documentation

    ## loading from cache

    ## see ?scRNAseq and browseVignettes('scRNAseq') for documentation

    ## loading from cache

    ## see ?scRNAseq and browseVignettes('scRNAseq') for documentation

    ## loading from cache

    ## see ?scRNAseq and browseVignettes('scRNAseq') for documentation

    ## loading from cache

    ## see ?scRNAseq and browseVignettes('scRNAseq') for documentation

    ## loading from cache

    ## see ?scRNAseq and browseVignettes('scRNAseq') for documentation

    ## loading from cache

    ## see ?scRNAseq and browseVignettes('scRNAseq') for documentation

    ## loading from cache

    ## see ?scRNAseq and browseVignettes('scRNAseq') for documentation

    ## loading from cache

    ## see ?scRNAseq and browseVignettes('scRNAseq') for documentation

    ## loading from cache

    ## see ?scRNAseq and browseVignettes('scRNAseq') for documentation

    ## loading from cache

    ## see ?scRNAseq and browseVignettes('scRNAseq') for documentation

    ## loading from cache

    ## see ?scRNAseq and browseVignettes('scRNAseq') for documentation

    ## loading from cache

    ## see ?scRNAseq and browseVignettes('scRNAseq') for documentation

    ## loading from cache

    ## see ?scRNAseq and browseVignettes('scRNAseq') for documentation

    ## loading from cache

    ## see ?scRNAseq and browseVignettes('scRNAseq') for documentation

    ## loading from cache

    ## see ?scRNAseq and browseVignettes('scRNAseq') for documentation

    ## loading from cache

    ## see ?scRNAseq and browseVignettes('scRNAseq') for documentation

    ## loading from cache

    ## see ?scRNAseq and browseVignettes('scRNAseq') for documentation

    ## loading from cache

Execute Targets workflow and load results

``` r
# Running the pipeline
preprocessed_seurat = run_targets_pipeline(
    input_data = input_data_path, 
    tissue = "pbmc",
    filter_input = TRUE, 
    RNA_assay_name = "originalexp", 
    sample_column = "Tissue"
) 
```

    ## Warning: Targets and globals must have unique names. Ignoring global objects
    ## that conflict with target names: sample_column, tissue, filter_input. Warnings
    ## like this one are important, but if you must suppress them, you can do so with
    ## Sys.setenv(TAR_WARN = "false").

    ## ▶ start target reference_file

    ## ● built target reference_file [0.388 seconds]

    ## ▶ start target reference_read

    ## ● built target reference_read [0.001 seconds]

    ## ▶ start target tissue_file

    ## ● built target tissue_file [0.002 seconds]

    ## ▶ start target tissue

    ## ● built target tissue [0 seconds]

    ## ▶ start target reference_label_fine

    ## ● built target reference_label_fine [0.002 seconds]

    ## ▶ start target reference_label_coarse

    ## ● built target reference_label_coarse [0 seconds]

    ## ▶ start target file

    ## ● built target file [0 seconds]

    ## ▶ start target filtered_file
    ## ● built target filtered_file [0.001 seconds]

    ## ▶ start target filter_input

    ## ● built target filter_input [0.001 seconds]

    ## ▶ start target read_file

    ## ● built target read_file [0.001 seconds]

    ## ▶ start branch input_read_e96ad026

    ## ● built branch input_read_e96ad026 [0.036 seconds]

    ## ● built pattern input_read

    ## ▶ start branch input_read_RNA_assay_762bd7f9

    ## ● built branch input_read_RNA_assay_762bd7f9 [0.031 seconds]

    ## ● built pattern input_read_RNA_assay

    ## ▶ start branch empty_droplets_tbl_585f9024
    ## ● built branch empty_droplets_tbl_585f9024 [8.776 seconds]

    ## ● built pattern empty_droplets_tbl

    ## ▶ start branch cell_cycle_score_tbl_7c4e4d58

    ## ● built branch cell_cycle_score_tbl_7c4e4d58 [0.219 seconds]

    ## ● built pattern cell_cycle_score_tbl

    ## ▶ start target sample_column_file

    ## ● built target sample_column_file [0.001 seconds]

    ## ▶ start target sample_column

    ## ● built target sample_column [0 seconds]

    ## ▶ start branch annotation_label_transfer_tbl_7c4e4d58

    ## ● built branch annotation_label_transfer_tbl_7c4e4d58 [20.34 seconds]

    ## ● built pattern annotation_label_transfer_tbl

    ## ▶ start branch alive_identification_tbl_18b1c13d

    ## ● built branch alive_identification_tbl_18b1c13d [1.52 seconds]

    ## ● built pattern alive_identification_tbl

    ## ▶ start branch non_batch_variation_removal_S_2a24bd5d

    ## ● built branch non_batch_variation_removal_S_2a24bd5d [5.653 seconds]

    ## ● built pattern non_batch_variation_removal_S

    ## ▶ start branch doublet_identification_tbl_2c729ad2

    ## ● built branch doublet_identification_tbl_2c729ad2 [11.837 seconds]

    ## ● built pattern doublet_identification_tbl

    ## ▶ start branch preprocessing_output_S_4fd096e6

    ## ● built branch preprocessing_output_S_4fd096e6 [0.279 seconds]

    ## ● built pattern preprocessing_output_S

    ## ▶ start target pseudobulk_preprocessing_SE

    ## ● built target pseudobulk_preprocessing_SE [5.426 seconds]

    ## ▶ end pipeline [1.366 minutes]

    ## Warning: 4 targets produced warnings. Run targets::tar_meta(fields = warnings,
    ## complete_only = TRUE) for the messages.

    ## HPCell says: you can read your output executing tar_read(preprocessing_output_S, store = "./")

``` r
# Load results
preprocessed_seurat
```

    ## $preprocessing_output_S_4fd096e6
    ## # A Seurat-tibble abstraction: 284 × 44
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
    ## # ℹ 38 more variables: nFeature_RNA <int>, percent.mito <dbl>,
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
