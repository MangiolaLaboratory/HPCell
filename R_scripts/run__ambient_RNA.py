#!/usr/bin/env python

import subprocess

# Install CellBender, pytorch and dependencies and create an environment for it according to instructions https://cellbender.readthedocs.io/en/latest/installation/index.html
# Activate conda environment and run CellBender command
# Make sure to start the script or submit the job with a GPU, otherwise pytorch won't work
# You can start the script on CPU only, for that comment out the --cuda, it will take much longer though
subprocess.run([
    "conda", "run", "-n", "CellBender",
    "cellbender", "remove-background",
    "--input", "/cellranger_output_folder/raw_feature_bc_matrix.h5",
    "--output", "/path_of_interest/cellbender_out_test.h5",
    "--cuda",
    "--expected-cells", "26000",
    "--total-droplets-included", "35000",
    "--fpr", "0.01",
    "--low-count-threshold", "150",
    "--epochs", "200"
], check=True)

#example to test
#subprocess.run([
#    "conda", "run", "-n", "CellBender",
#    "cellbender", "remove-background",
#    "--input", "/stornext/Bioinf/data/bioinf-data/Papenfuss_lab/projects/mangiola.s/PostDoc/covid19pbmc/data/all_batches/raw_counts/C120_COVID_PBMC_batch3/extdata/CellRanger/C120_batch3_1/outs/multi/count/raw_feature_bc_matrix.h5",
#    "--output", "/stornext/Bioinf/data/bioinf-data/Papenfuss_lab/projects/mangiola.s/PostDoc/covid19pbmc/data/all_batches/raw_counts/C120_COVID_PBMC_batch3/batch3_1_cellbender_out_test.h5",
#    "--cuda",
#    "--expected-cells", "26000",
#    "--total-droplets-included", "35000",
#    "--fpr", "0.01",
#    "--low-count-threshold", "150",
#    "--epochs", "200"
#], check=True)
