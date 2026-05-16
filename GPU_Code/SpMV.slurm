#!/bin/bash
#SBATCH --partition=edu-short
#SBATCH --account=gpu.computing26
#SBATCH --time=00:4:30
#SBATCH --gres=gpu:1
#SBATCH --tasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8

#SBATCH --job-name=SpMV
#SBATCH --output=program_output/R-%x.%j.out
#SBATCH --error=program_output/R-%x.%j.err

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

# 3. Optional: Bind threads to cores for better performance 
export OMP_PLACES=cores
export OMP_PROC_BIND=close

module load CUDA/12.5.0

for MATRIX_NUMBER in $(seq 1 11); do
    ./bin/SpMV "$MATRIX_NUMBER"
done
