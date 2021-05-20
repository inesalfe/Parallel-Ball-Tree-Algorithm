#!/bin/bash

#SBATCH --job-name=ballAlgIR
#SBATCH --output=pts2.txt
#SBATCH --error=stderr2
#SBATCH --ntasks=4
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
srun ballAlg 3 5000000 0