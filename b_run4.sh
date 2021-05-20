#!/bin/bash

#SBATCH --job-name=ballAlgIR
#SBATCH --output=pts4.txt
#SBATCH --error=stderr4
#SBATCH --ntasks=4
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
srun ballAlg 3 20000000 0