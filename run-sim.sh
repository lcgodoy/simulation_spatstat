#!/bin/bash
#SBATCH --partition=general       # you might have to change the partition
#SBATCH --job-name=spatstat_paper
#SBATCH --array=1-1000%20         # run it a thousand times "20 by 20"
#SBATCH --ntasks=1

module purge
module load singularity/3.5.2

## avoiding implicit paralellism
export OMP_NUM_THREADS=1

singularity exec spatstat.sif \
	    Rscript --vanilla scripts/01-sim-data_parallel.R \
	    $SLURM_ARRAY_TASK_ID

singularity exec spatstat.sif \
	    Rscript --vanilla scripts/02-ht_parallel.R \
	    $SLURM_ARRAY_TASK_ID
