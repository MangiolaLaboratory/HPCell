# HPCell

## Analysis automated pipeline

Flowchart
https://app.mural.co/t/covid7029/m/covid7029/1656652076667/c47e104697d76b36b8dee3bd05d8de9d96d99efd?sender=udbe1ea99c9618abb07196978

Clone the repository

```{bash}

git clone git@github.com:susansjy22/jascap.git

```

Enter in the Jascap directory within R and activate renv for reproducibility

```{bash}

module load R/4.2.1

```

Install Jascap package 

```{bash}

remote::install_github("git@github.com:susansjy22/jascap.git")

```

load jascap package 

```{bash}

library(jascap)

```

load input and reference data

```{bash}

# Load input data (can be a list of directories or single directory)
input_data <- c("/home/users/allstaff/si.j/test_jascap/input/CB150T04X__batch14.rds","/home/users/allstaff/si.j/test_jascap/input/CB291T01X__batch8.rds")

#Load reference data 
input_reference <- "/home/users/allstaff/si.j/jascap/data/Data/jiayi_files/reference_azimuth.rds"

```

Execute Targets workflow and load results

```{bash}
#Create store directory 
store =  tempfile(tmpdir = "/stornext/General/scratch/GP_Transfer/si.j")

#Execute pipeline

preprocessed_seurat = run_targets_pipeline(input_data, store , input_reference)

#Load results

preprocessed_seurat

```



