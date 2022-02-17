#!/bin/bash
#SBATCH --partition=general                # Name of Partition
#SBATCH --job-name=spatstat_paper
#SBATCH --array=1000%20                    # run it a thousand times "20 by 20"
#SBATCH --ntasks=1
#SBATCH --mail-type=END,FAIL               # Event(s) that triggers email notification (BEGIN,END,FAIL,ALL)

module purge
module load singularity/3.5.2

## avoiding implicit paralellism
export OMP_NUM_THREADS=1

singularity exec spatstat.sif \
	    Rscript --vanilla scripts/01-sim-data.R \
	    $SLURM_ARRAY_TASK_ID

singularity exec spatstat.sif \
	    Rscript --vanilla scripts/02-ht.R \
	    $SLURM_ARRAY_TASK_ID
