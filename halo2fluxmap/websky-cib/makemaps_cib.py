#!/bin/bash -l
#SBATCH -q premium
#SBATCH -N 16
#SBATCH --tasks-per-node=8
#SBATCH -C haswell
#SBATCH -A mp107
#SBATCH -t 8:00:00
#SBATCH -J h2fm
#SBATCH -o h2fm.o%j
#SBATCH -L SCRATCH,project

export OMP_PROC_BIND=true
export OMP_PLACES=threads
export OMP_NUM_THREADS=8

# run driver
conda activate halo2fluxmap
srun -n $SLURM_NTASKS --mpi=pmi2 python scripts/makemaps_cib_viero.py
