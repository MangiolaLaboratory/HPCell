# HPCell

## Analysis automated pipeline

Flowchart
https://app.mural.co/t/covid7029/m/covid7029/1656652076667/c47e104697d76b36b8dee3bd05d8de9d96d99efd?sender=udbe1ea99c9618abb07196978

Clone the repository

```{r}

git clone git@github.com:susansjy22/HPCell.git

```

Install Jascap package 

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
file_path = tempfile(tmpdir = "~/Documents") |> paste0(".rds")
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
store =  tempfile(tmpdir = "~/Documents")

#Execute pipeline
#Tissue modalities: pbmc, solid, atypical
preprocessed_seurat = run_targets_pipeline(c(file_path,file_path), store ,input_reference_path, tissue = "pbmc")

#Load results

preprocessed_seurat

```



