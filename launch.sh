#!/bin/bash

#SBATCH --job-name=A20-k1-ginf-l3
#SBATCH --partition=short
#SBATCH --time=02:00:00
#SBATCH --nodes=1           ### Number of Nodes
#SBATCH --ntasks-per-node=1 ### Number of tasks (MPI processes)
#SBATCH --cpus-per-task=128
#SBATCH --mem=16G
#SBATCH --dependency=afterany:5932351
#SBATCH --signal=15@180
#SBATCH --output=job_%j.out
#SBATCH --error=job_%j.err

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OMP_PLACES=threads
export OMP_BIND_PROC=true

# ml intel-modules intel intel.mpi
ml amd-modules

# Start 'myprog' with input from bucket, 
# and output to our temporary directory 
srun ball0x.exe > log.out 2> err.out 

