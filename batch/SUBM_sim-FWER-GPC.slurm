#!/bin/bash

#SBATCH --job-name=sim-FWER-GPC  # Job name
#SBATCH --array=1-200     # Number of jobs
#SBATCH --partition=standard    # Queue. Use long for run time >8h and standard otherwise
#SBATCH --mem=1G # Memory limit, e.g. reserve 1 GB memory 
#SBATCH --output=output/sim-FWER-GPC/slurm-%x-%a-%J.out # Direct output to subdirectory
#SBATCH --error=output/sim-FWER-GPC/slurm-%x-%a-%J.out # Same file for output and error

R CMD BATCH --vanilla BATCH_sim-FWER-GPC.R output/sim-FWER-GPC/R-$SLURM_JOB_NAME-$SLURM_ARRAY_TASK_ID-$SLURM_JOB_ID.Rout

## go to directory    ## cd x:/GPC/book/
                      ## cd /projects/biostat01/people/hpl802/GPC/book/
## clean outputs      ## rm -r ./output/sim-FWER-GPC/*
## clean results      ## rm -r ./Results/sim-FWER-GPC/*
## submission command ## sbatch SUBM_sim-FWER-GPC.slurm

## submission output  ## Submitted batch job 4371973 (time using C-c t: onsdag 23-09-27 at 17:41) 

## init               ## module module load gcc/11.2.0 R/4.1.2
                      ## cat .bash_profile (in the root directory on the server to visualize the current configuration)
                      ## vi .bash_profile (in the root directory on the server to modify the current configuration)

## commands           ## squeue           : view current jobs 
                      ## squeue -u id     : view current jobs for user id
                      ## scancel          : delete the job with job id 1034 from the queue type
                      ## sinfo            : view state of the queues
                      ## sshare           : check own usage of the server
                        
## documentation      ## from SLURM: https://slurm.schedmd.com/pdfs/summary.pdf: key functions
                      ##           : https://slurm.schedmd.com/sbatch.html (filename pattern)
                      ## from KU   : https://hpc.ku.dk/documentation/slurm.html
                      ##           : https://github.com/aejensen/computing-biostat-ucph/blob/main/README.md
