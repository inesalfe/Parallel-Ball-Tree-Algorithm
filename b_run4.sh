#!/bin/bash

#SBATCH --job-name=ballAlg_i_r
#SBATCH --output=pts4.txt
#SBATCH --error=stderr4
#SBATCH --ntasks=16
#SBATCH --ntasks-per-node=1
srun ballAlg 3 20000000 0