#!/bin/bash

#SBATCH --job-name=ballAlgIR
#SBATCH --error=stderr.txt
#SBATCH --output=pts.txt
#SBATCH --ntasks=4
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OMP_NESTED=true

loop_variable=10

# 20 1000000 0

echo "--- 20 200 0 ---"

sum=0
for (( i=0; i<loop_variable; i++ ))
do
	time=`srun ballAlg 20 200 0 2>&1`
	wait
	echo "$time"
	sum=$(echo $sum + $time | bc -l);
done

avg=$(echo $sum / $loop_variable | bc -l);
avg=`printf "%.1f" $avg`
echo "Average time: $avg"