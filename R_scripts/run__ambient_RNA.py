#!/bin/bash
#SBATCH --job-name=cellbender 
#SBATCH --time=48:00:00
#SBATCH --ntasks=1
#SBATCH --mem=7G
#SBATCH --partition=gpuq
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=2
#SBATCH -o slurm.%j.out
#SBATCH -e slurm.%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=odainic.a@wehi.edu.au
module load anaconda3/latest
conda activate /home/users/allstaff/odainic.a/.conda/envs/CellBender
cellbender remove-background \
     --input /stornext/Bioinf/data/bioinf-data/Papenfuss_lab/projects/mangiola.s/PostDoc/covid19pbmc/data/all_batches/raw_counts/C120_COVID_PBMC_batch3/extdata/CellRanger/C120_batch3_1/outs/multi/count/raw_feature_bc_matrix.h5 \
                 --output /stornext/Bioinf/data/bioinf-data/Papenfuss_lab/projects/odainic.a/CellBender/mytest/batch3_1_cellbender_out.h5 \
                 --cuda \
                 --expected-cells 26000 \
                 --total-droplets-included 40000 \
                 --fpr 0.01 \
                 --low-count-threshold 150 \
                 --epochs 200
conda deactivate
