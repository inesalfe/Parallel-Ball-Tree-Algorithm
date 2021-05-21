#!/bin/bash

loop_variable=10
proc="8"

echo "--- 20 1000000 0 ---"

sum=0
for (( i=0; i<loop_variable; i++ ))
do
	time=`srun --ntasks-per-node=1 -n "$proc" ballAlg 20 1000000 0 2>&1`
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
	time=`srun --ntasks-per-node=1 -n "$proc" ballAlg 3 5000000 0 2>&1`
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
	time=`srun --ntasks-per-node=1 -n "$proc" ballAlg 4 10000000 0 2>&1`
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
	time=`srun --ntasks-per-node=1 -n "$proc" ballAlg 3 20000000 0 2>&1`
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
	time=`srun --ntasks-per-node=1 -n "$proc" ballAlg 4 20000000 0 2>&1`
	wait
	echo "$time"
	sum=`echo "$sum + $time" | bc -l`
done

avg=`echo "$sum / $loop_variable" | bc -l`
avg=`printf "%.1f" $avg`
echo "Average time: $avg"