#!/bin/bash

#---------------------------------------------------------------------------------
# Account information

#SBATCH --account=pi-bata0              # basic (default), staff, phd, faculty
#SBATCH --qos=glass
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=Yuwei.Zhou@chicagobooth.edu

#---------------------------------------------------------------------------------
# Resources requested

#SBATCH --partition=standard       # standard (default), long, gpu, mpi, highmem
#SBATCH --cpus-per-task=1          # number of CPUs requested (for parallel tasks)
#SBATCH --mem=2G           # requested memory
#SBATCH --time=0-72:00:00          # wall clock limit (d-hh:mm:ss)
#SBATCH --ntasks=1

#---------------------------------------------------------------------------------
# Job specific name (helps organize and track progress of jobs)

#SBATCH --job-name=my_batch_job    # user-defined job name

#---------------------------------------------------------------------------------
# Print some useful variables

echo "Job ID: $SLURM_JOB_ID"
echo "Job User: $SLURM_JOB_USER"
echo "Num Cores: $SLURM_JOB_CPUS_PER_NODE"

#---------------------------------------------------------------------------------
# Load necessary modules for the job

date
pwd

g++ BerRF.cpp -o BF
srun BF

#---------------------------------------------------------------------------------
# Commands to execute below...
