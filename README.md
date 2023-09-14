# HPCell

## Analysis automated pipeline

Flowchart

https://app.mural.co/t/covid7029/m/covid7029/1656652076667/c47e104697d76b36b8dee3bd05d8de9d96d99efd?sender=udbe1ea99c9618abb07196978

Clone the repository

```{r}

git clone git@github.com:susansjy22/HPCell.git

```

Install Jascap/HPCell package 

```{r}

remote::install_github("git@github.com:susansjy22/HPCell.git")

```

load jascap package 

```{r}

library(jascap)

```

load input and reference data


```{r}

# Load input data (can be a list of directories or single directory)

library(Seurat)
library(scRNAseq)

single_cell_data = ChenBrainData(ensembl=FALSE,location=FALSE)

file_path = tempfile(tmpdir = "~/HPCell") |> paste0(".rds")

# Convert from SingleCellExperiment to Seurat class
single_cell_data |> Seurat::as.Seurat(data = NULL, assay = NULL) |> saveRDS(file_path)
input_data_seurat <- readRDS(file_path)

# # Get current default assay
# DefaultAssay(object = input_data_seurat)
# 
# # If not RNA 
# original_assay<- DefaultAssay(object = input_data_seurat)
# 
# # Create RNA assay to demo switching default assays
# new.assay <- input_data_seurat[[original_assay]]
# Key(object = new.assay) <- "RNA_"
# input_data_seurat[["RNA"]] <- new.assay
# # switch default assay to RNA2
# DefaultAssay(object = input_data_seurat) <- "RNA"
# DefaultAssay(object = input_data_seurat)

#Load reference data 

library(Azimuth)

library(SeuratData)

InstallData("pbmcsca")

input_reference_path = "~/HPCell/pbmcsca.rds"
pbmcsca.updated = UpdateSeuratObject(object = pbmcsca)
pbmcsca.updated |> saveRDS(input_reference_path)

```

Execute Targets workflow and load results


```{r}

#Create store directory 

store =  tempfile(tmpdir = "~/stornext/General/scratch/GP_Transfer/si.j/test_sc_postmod")

store <- "~/stornext/General/scratch/GP_Transfer/si.j/test_sc_postmod"


#Execute pipeline

#Tissue types: pbmc, solid, atypical

preprocessed_seurat = run_targets_pipeline(
    input_data = file_path, 
    store =  "~/stornext/General/scratch/GP_Transfer/si.j/test_single_cell_data", 
    input_reference = input_reference_path,
    tissue= "pbmc",
    computing_resources = crew_controller_local(workers = 1), 
    debug_step = "empty_droplets_tbl"
  )


#Load results

preprocessed_seurat

```

Debugging ~R/functions.R

If an error occurs in one of the Target objects whilst executing the pipeline, you can debug by going inside the function to see where the error occurred

Step-by-step guide

1. Set the debug_step argument to the name of the failed Targets object 
   - Example: "annotation_label_transfer_tbl"

```{r}

run_targets_pipeline(
    input_data = file_path, 
    store =  tempfile(tmpdir = "."), 
    input_reference,
    tissue,
    computing_resources = crew_controller_local(workers = 1), 
    debug_step = "annotation_label_transfer_tbl"
  )

```

2. Re-run the pipeline; Pipeline should halt execution at the specified target and enter into debug mode 

3. Use debugonce() command and specify the function matching the Targets object 

```{r}

debugonce(annotation_label_transfer)

```
4. Locate the function and associated code within ~R/functions.R Rscript

5. Debug as normal 

!!!IMPORTANT!!! 
   
   Remember to replace the function arguments with arguments provided in ~R/execute_pipeline.R file

