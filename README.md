# PBMC

## Analysis automated pipeline

Clone the repository

```{bash}

git clone git@github.com:Melbourne-COVID-Predict/jascap.git

```

Enter in the jascap directory within R and activate renv, for reproducibility

```{r}
install.packages("renv")
 renv::activate()
 renv::restore()
```

Install makeflow

```{bash}
# Load miniconda
module load miniconda3

# run this command once
$ conda create -n cctools-env -y -c conda-forge --strict-channel-priority python ndcctools

# run this command every time you want to use cctools
$ conda activate cctools-env

```

Create input and result directories 

```{bash}

# Define dataset directory, for example test_pipeline
master_directory=test_jacap
mkdir $master_directory
result_directory=$master_directory/results
mkdir $result_directory
input_directory=$master_directory/input
mkdir $input_directory

# The R directory in the same location where you cloned jascap
code_directory=~/PostDoc/jascap

# This is the location of the metadata, which should match by `sample` with the column in the seurat objects
metadata_path=~/PostDoc/covid19pbmc/data/3_prime_batch_1/metadata.rds

# This is in a shared location 
reference_azimuth_path=/stornext/Bioinf/data/bioinf-data/Papenfuss_lab/projects/reference_azimuth.rds
```
Add input files in the input directory. 

1) Each input file should be a seurat object
2) Each input file should just inbclude row counts in the `RNA` assay
3) Each input file should include only one sample per file
4) Each input file should be named `SAMPLE_NAME`.rds
5) All `SAMPLE_NAME` should be unique for each sample (if a sample has different timepoints the `SAMPLE_NAME` should be made unique, for example `SAMPLE_NAME_TIME_POINT`)

You should copy those files in your input directory, for example 

```{bash}
cp myfiles/* $input_directory
```

Execute pipeline using `makeflow`

```{bash}
# Create $input_directory/pipeline.makeflow
Rscript $code_directory/R/create_pipeline_makefile.R $result_directory $input_directory $code_directory $metadata_path $reference_azimuth_path

# Execute makeflow
conda activate cctools-env
makeflow -J 200 -T slurm $result_directory/pipeline.makeflow 
```

Monitor progress

```{bash}
makeflow_monitor $result_directory/pipeline.makeflow.makeflowlog
```
