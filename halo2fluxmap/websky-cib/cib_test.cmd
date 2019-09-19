#!/bin/bash
# parallel job using 16 processors. and runs for 4 hours (max)
#SBATCH --nodes=1 # node count
#SBATCH --ntasks-per-node=20
#SBATCH --mem 680G
#SBATCH --time=1:00:00 # sends mail when process begins, and
# when it ends. Make sure you define your email
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH --mail-user=zq@princeton.edu
#
# Load openmpi environment
#
module load anaconda intel intel-mkl intel-mpi mpi4py
cd /home/zequnl/xgpaint/halo2fluxmap/websky-cib
srun -n $SLURM_NTASKS makemaps_cib.py 

