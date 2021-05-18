#!/bin/sh

loop_variable=10
os="0"
proc="2"

# 20 1000000 0

echo "--- 20 1000000 0 ---"

sum=0
for (( i=0; i<loop_variable; i++ ))
do
	time=`srun -n "$proc" ballAlg 20 1000000 "$os" 2>&1`
	echo "$time"
	sum=$(echo $sum + $time | bc -l);
done

avg=$(echo $sum / $loop_variable | bc -l);
avg=`printf "%.1f" $avg`
echo "Average time: $avg"

# 3 5000000 0

echo "--- 3 5000000 0 ---"

sum=0
for (( i=0; i<loop_variable; i++ ))
do
	time=`srun -n "$proc" ballAlg 3 5000000 "$os" 2>&1`
	echo "$time"
	sum=$(echo $sum + $time | bc -l);
done

avg=$(echo $sum / $loop_variable | bc -l);
avg=`printf "%.1f" $avg`
echo "Average time: $avg"

# 4 10000000 0

echo "--- 4 10000000 0 ---"

sum=0
for (( i=0; i<loop_variable; i++ ))
do
	time=`srun -n "$proc" ballAlg 4 10000000 "$os" 2>&1`
	echo "$time"
	sum=$(echo $sum + $time | bc -l);
done

avg=$(echo $sum / $loop_variable | bc -l);
avg=`printf "%.1f" $avg`
echo "Average time: $avg"

# 3 20000000 0

echo "--- 3 20000000 0 ---"

sum=0
for (( i=0; i<loop_variable; i++ ))
do
	time=`srun -n "$proc" ballAlg 3 20000000 "$os" 2>&1`
	echo "$time"
	sum=$(echo $sum + $time | bc -l);
done

avg=$(echo $sum / $loop_variable | bc -l);
avg=`printf "%.1f" $avg`
echo "Average time: $avg"

# 4 20000000 0

echo "--- 4 20000000 0 ---"

sum=0
for (( i=0; i<loop_variable; i++ ))
do
	time=`srun -n "$proc" ballAlg 4 20000000 "$os" 2>&1`
	echo "$time"
	sum=$(echo $sum + $time | bc -l);
done

avg=$(echo $sum / $loop_variable | bc -l);
avg=`printf "%.1f" $avg`
echo "Average time: $avg"