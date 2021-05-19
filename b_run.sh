#!/bin/bash

#SBATCH --job-name=ballAlg_i_r
#SBATCH --output=pts.txt
#SBATCH --error=stderr
#SBATCH --ntasks=4
#SBATCH --ntasks-per-node=1
#SBATCH -–cpus-per-task=4
srun ballAlg 3 50 0