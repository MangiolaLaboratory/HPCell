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

single_cell_data |> Seurat::as.Seurat(data = NULL) |> saveRDS(file_path)


#Load reference data 

library(Azimuth)

library(SeuratData)

InstallData("pbmcsca")

input_reference_path = "reference_azimuth.rds"

LoadData("pbmcsca") |> saveRDS(input_reference_path)

```

Execute Targets workflow and load results


```{r}

#Create store directory 

store =  tempfile(tmpdir = "~/stornext/General/scratch/GP_Transfer/si.j/test_single_cell_data")

store <- "~/stornext/General/scratch/GP_Transfer/si.j/test_single_cell_data"


#Execute pipeline

#Tissue types: pbmc, solid, atypical

preprocessed_seurat = run_targets_pipeline(input_data = "/home/users/allstaff/si.j/HPCell/file29c24715b278.rds", store ,input_reference_path, tissue = "pbmc")


#Load results

preprocessed_seurat

```

Debugging ~R/functions.R

If an error occurs in one of the Target objects whilst executing the pipeline, you can debug by going inside the function to see where the error occurred

Step-by-step guide

1. Set the debug_step argument to the name of the failed Targets object (e.g., "annotation_label_transfer_tbl")

2. Re-run the pipeline; Pipeline should halt execution at the specified target and enter into debug mode 

3. Use the debugonce() command for function matching the Targets object 

```{r}

debugonce(annotation_label_transfer)

```
4. Locate the function and associated code within ~R/functions.R Rscript

5. Debug as normal 

!!!IMPORTANT!!! 
   
   Remember to replace the function arguments with arguments provided in ~R/execute_pipeline.R file

