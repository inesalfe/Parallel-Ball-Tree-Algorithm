#!/bin/bash

#SBATCH --job-name=ballAlg_i_r
#SBATCH --output=pts5.txt
#SBATCH --error=stderr5
#SBATCH --ntasks=16
#SBATCH --ntasks-per-node=1
srun ballAlg 4 20000000 0