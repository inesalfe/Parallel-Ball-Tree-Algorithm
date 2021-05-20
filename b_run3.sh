#!/bin/bash

#SBATCH --job-name=ballAlgIR
#SBATCH --output=pts3.txt
#SBATCH --error=stderr3
#SBATCH --ntasks=4
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
srun ballAlg 4 10000000 0