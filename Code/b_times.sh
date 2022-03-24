#!/bin/bash

#SBATCH --job-name=ballAlgIR
#SBATCH --error=stderr_time.txt
#SBATCH --output=times.txt
#SBATCH --ntasks=4
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --exclude=lab6p[1-9]

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OMP_NESTED=true

loop_variable=10

# echo "--- 2 5 0 ---"

# sum=0
# for (( i=0; i<loop_variable; i++ ))
# do
# 	time=`srun ballAlg 2 5 0 2>&1`
# 	wait
# 	echo "$time"
# 	sum=`echo "$sum + $time" | bc -l`
# done

# avg=`echo "$sum / $loop_variable" | bc -l`
# avg=`printf "%.1f" $avg`
# echo "Average time: $avg"

# echo "--- 2 8 0 ---"

# sum=0
# for (( i=0; i<loop_variable; i++ ))
# do
# 	time=`srun ballAlg 2 8 0 2>&1`
# 	wait
# 	echo "$time"
# 	sum=`echo "$sum + $time" | bc -l`
# done

# avg=`echo "$sum / $loop_variable" | bc -l`
# avg=`printf "%.1f" $avg`
# echo "Average time: $avg"

# echo "--- 3 3 3 ---"

# sum=0
# for (( i=0; i<loop_variable; i++ ))
# do
# 	time=`srun ballAlg 3 3 3 2>&1`
# 	wait
# 	echo "$time"
# 	sum=`echo "$sum + $time" | bc -l`
# done

# avg=`echo "$sum / $loop_variable" | bc -l`
# avg=`printf "%.1f" $avg`
# echo "Average time: $avg"

# echo "--- 5 4 3 ---"

# sum=0
# for (( i=0; i<loop_variable; i++ ))
# do
# 	time=`srun ballAlg 5 4 3 2>&1`
# 	wait
# 	echo "$time"
# 	sum=`echo "$sum + $time" | bc -l`
# done

# avg=`echo "$sum / $loop_variable" | bc -l`
# avg=`printf "%.1f" $avg`
# echo "Average time: $avg"

# echo "--- 2 6 8 ---"

# sum=0
# for (( i=0; i<loop_variable; i++ ))
# do
# 	time=`srun ballAlg 2 6 8 2>&1`
# 	wait
# 	echo "$time"
# 	sum=`echo "$sum + $time" | bc -l`
# done

# avg=`echo "$sum / $loop_variable" | bc -l`
# avg=`printf "%.1f" $avg`
# echo "Average time: $avg"

echo "--- 20 1000000 0 ---"

sum=0
for (( i=0; i<loop_variable; i++ ))
do
	time=`srun ballAlg 20 1000000 0 2>&1`
	wait
	echo "$time"
	sum=`echo "$sum + $time" | bc -l`
done

avg=`echo "$sum / $loop_variable" | bc -l`
avg=`printf "%.1f" $avg`
echo "Average time: $avg"

echo "--- 3 5000000 0 ---"

sum=0
for (( i=0; i<loop_variable; i++ ))
do
	time=`srun ballAlg 3 5000000 0 2>&1`
	wait
	echo "$time"
	sum=`echo "$sum + $time" | bc -l`
done

avg=`echo "$sum / $loop_variable" | bc -l`
avg=`printf "%.1f" $avg`
echo "Average time: $avg"

echo "--- 4 10000000 0 ---"

sum=0
for (( i=0; i<loop_variable; i++ ))
do
	time=`srun ballAlg 4 10000000 0 2>&1`
	wait
	echo "$time"
	sum=`echo "$sum + $time" | bc -l`
done

avg=`echo "$sum / $loop_variable" | bc -l`
avg=`printf "%.1f" $avg`
echo "Average time: $avg"

echo "--- 3 20000000 0 ---"

sum=0
for (( i=0; i<loop_variable; i++ ))
do
	time=`srun ballAlg 3 20000000 0 2>&1`
	wait
	echo "$time"
	sum=`echo "$sum + $time" | bc -l`
done

avg=`echo "$sum / $loop_variable" | bc -l`
avg=`printf "%.1f" $avg`
echo "Average time: $avg"

echo "--- 4 20000000 0 ---"

sum=0
for (( i=0; i<loop_variable; i++ ))
do
	time=`srun ballAlg 4 20000000 0 2>&1`
	wait
	echo "$time"
	sum=`echo "$sum + $time" | bc -l`
done

avg=`echo "$sum / $loop_variable" | bc -l`
avg=`printf "%.1f" $avg`
echo "Average time: $avg"