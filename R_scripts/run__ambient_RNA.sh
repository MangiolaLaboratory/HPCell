:'
Install CellBender and create the environment before using this script
1a. If you are on a HPC, load anaconda with module load
1b. If you are on Linux and you do not have anaconda, install it following the steps https://docs.anaconda.com/free/anaconda/install/linux/
2. Create an anaconda environment where CellBender will be installed:
$ conda create -n CellBender python=3.7
3. Activate the environment to install the packages inside:
$ conda activate CellBender
4. Install the pytables module:
(CellBender) $ conda install -c anaconda pytables
5. Install pytorch to be able to use CUDA:
(CellBender) $ conda install pytorch torchvision torchaudio pytorch-cuda=11.7 -c pytorch -c nvidia
This might deviate based on what system you are, what GPU you have and what drivers you are using. The upper code is for Linux installation through conda and python language. For other pytorch versions check https://pytorch.org/get-started/locally/
6. Clone the CellBender repository from github and install CellBender:
(CellBender) $ git clone https://github.com/broadinstitute/CellBender.git
(CellBender) $ pip install -e CellBender
7. You are all set to start, for CellBender arguments and troubleshooting check https://cellbender.readthedocs.io/en/latest/help_and_reference/remove_background/index.html
When starting CellBender, you activate your conda environment with conda activate CellBender. Then you run cellbender remove-background --input --output and other arguments
If starting a slurm job, make sure to load anaconda and activate the environment for CellBender before starting. See the below example
'
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
                 --output /stornext/Bioinf/data/bioinf-data/Papenfuss_lab/projects/odainic.a/CellBender/mytest/batch3_1_cellbender_out_test.h5 \
                 --cuda \
                 --expected-cells 26000 \
                 --total-droplets-included 40000 \
                 --fpr 0.01 \
                 --low-count-threshold 150 \
                 --epochs 200
conda deactivate
