#!/bin/bash

#SBATCH --job-name=ballAlg_i_r
#SBATCH --output=pts2.txt
#SBATCH --error=stderr2
#SBATCH --ntasks=32
#SBATCH --ntasks-per-node=1
srun ballAlg 3 5000000 0