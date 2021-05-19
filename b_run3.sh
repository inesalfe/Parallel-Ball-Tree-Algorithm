#!/bin/bash

#SBATCH --job-name=ballAlg_i_r
#SBATCH --output=pts3.txt
#SBATCH --error=stderr3
#SBATCH --ntasks=32
#SBATCH --ntasks-per-node=1
srun ballAlg 4 10000000 0