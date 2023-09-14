# HPCell

## Analysis automated pipeline

Flowchart
https://app.mural.co/t/covid7029/m/covid7029/1656652076667/c47e104697d76b36b8dee3bd05d8de9d96d99efd?sender=udbe1ea99c9618abb07196978

Clone the repository

```{bash}

git clone git@github.com:susansjy22/HPCell.git

```

Install Jascap package 

```{bash}

remote::install_github("git@github.com:susansjy22/HPCell.git")

```

load jascap package 

```{bash}

library(jascap)

```

load input and reference data

```{bash}

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

```{bash}
#Create store directory 
store =  tempfile(tmpdir = "/stornext/General/scratch/GP_Transfer/si.j")
store <- "~/stornext/General/scratch/GP_Transfer/si.j/test_single_cell_data"
#Execute pipeline

preprocessed_seurat = run_targets_pipeline(input_data, store , input_reference)

#Load results

preprocessed_seurat

```



