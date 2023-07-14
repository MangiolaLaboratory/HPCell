# PBMC

## Analysis automated pipeline

Flowchart
https://app.mural.co/t/covid7029/m/covid7029/1656652076667/c47e104697d76b36b8dee3bd05d8de9d96d99efd?sender=udbe1ea99c9618abb07196978

Clone the repository

```{bash}

git clone git@github.com:Melbourne-COVID-Predict/jascap.git

```

Enter in the jascap directory within R and activate renv, for reproducibility

```{r}
module load R/4.2.1
```
open R in the terminal:

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
reports_directory=$master_directory/results/reports
mkdir $result_directory
input_directory=$master_directory/input
mkdir $input_directory

# The R directory in the same location where you cloned jascap
code_directory=~/PostDoc/jascap

# This is the location of the metadata, which should match by `sample` with the column in the seurat objects. The metadata file needs to be in a different directory than the input file.
metadata_path=~/PostDoc/covid19pbmc/data/3_prime_batch_1/metadata.rds

# This is in a shared location 
reference_azimuth_path=/stornext/Bioinf/data/bioinf-data/Papenfuss_lab/projects/reference_azimuth.rds
```

>*Note: If no `metadata.rds` is available: create a data frame with 2 columns (sample | batch) with the `samples` matching the `SAMPLE_NAME` in the input file, >as described below. Set `batch` as Sample_name.*

| sample | batch |
| :---: | :---: |
| spleen | spleen|
| liver | liver |

Add input files in the input directory. 

1) Each input file should be a seurat object
2) Each input file should just include row counts in the `RNA` assay
3) Each input file should include only one sample 
4) Each input file should be named `SAMPLE_NAME`.rds
5) All `SAMPLE_NAME` should be unique for each sample (if a sample has different time points the `SAMPLE_NAME` should be made unique, for example `SAMPLE_NAME_TIME_POINT`)

You should copy those files in your input directory, for example 

```{bash}
cp myfiles/* $input_directory
```

Execute pipeline using `makeflow`

*The available modalities are*
- preprocessing (filtering and annotation without integration)
- fast (preprocessing + differential transcription, tissue composition, and cell communication analyses)
- complete (still not ready)

*The available tissue annotation are*
- pbmc (preprocessing based on Seurat Azimuth pbmc annotation)
- solid (preprocessing based on SingleR blueprint annotation)
- atypical (preprocessing based on unsupervised clustering)

```{bash}
# Create $input_directory/pipeline.makeflow
#
# The pipeline has 4 modes (first argument of create_pipeline_makefile.R)
# preprocessing, fast_pipeline, slow_pipeline, complete
Rscript $code_directory/R_scripts/create_pipeline_makefile.R complete pbmc unfiltered $result_directory $reports_directory $input_directory $code_directory $metadata_path $reference_azimuth_path

# Execute makeflow
conda activate cctools-env
makeflow -J 200 -T slurm $result_directory/pipeline.makeflow 
```

Monitor progress

```{bash}
makeflow_monitor $result_directory/pipeline.makeflow.makeflowlog
```
