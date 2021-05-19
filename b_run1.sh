#!/bin/bash

#SBATCH --job-name=ballAlg_i_r
#SBATCH --output=pts1.txt
#SBATCH --error=stderr1
#SBATCH --ntasks=32
#SBATCH --ntasks-per-node=1
srun ballAlg 20 1000000 0