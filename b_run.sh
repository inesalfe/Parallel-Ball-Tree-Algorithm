#!/bin/bash

#SBATCH --job-name=ballAlgIR
#SBATCH --output=pts.txt
#SBATCH --error=stderr
#SBATCH --ntasks=4
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
srun ballAlg 20 1000000 0