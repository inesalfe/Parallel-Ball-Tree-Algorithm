#!/bin/sh

loop_variable=10
os="1"

# 2 5 1

echo "--- 2 5 0 ---"

sum=0
for (( i=0; i<loop_variable; i++ ))
do
	time=`./ballAlg 2 5 "$os" 2>&1`
	echo "$time"
	sum=$(echo $sum + $time | bc -l);
done

avg=$(echo $sum / $loop_variable | bc -l);
avg=`printf "%.1f" $avg`
echo "Average time: $avg"

# 2 8 1

echo "--- 2 8 0 ---"

sum=0
for (( i=0; i<loop_variable; i++ ))
do
	time=`./ballAlg 2 8 "$os" 2>&1`
	echo "$time"
	sum=$(echo $sum + $time | bc -l);
done

avg=$(echo $sum / $loop_variable | bc -l);
avg=`printf "%.1f" $avg`
echo "Average time: $avg"

# 20 1000000 0

echo "--- 20 1000000 0 ---"

sum=0
for (( i=0; i<loop_variable; i++ ))
do
	time=`./ballAlg 20 1000000 "$os" 2>&1`
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
	time=`./ballAlg 3 5000000 "$os" 2>&1`
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
	time=`./ballAlg 4 10000000 "$os" 2>&1`
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
	time=`./ballAlg 3 20000000 "$os" 2>&1`
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
	time=`./ballAlg 4 20000000 "$os" 2>&1`
	echo "$time"
	sum=$(echo $sum + $time | bc -l);
done

avg=$(echo $sum / $loop_variable | bc -l);
avg=`printf "%.1f" $avg`
echo "Average time: $avg"

# 3 3 3

echo "--- 3 3 3 ---"

sum=0
for (( i=0; i<loop_variable; i++ ))
do
	time=`./ballAlg 3 3 3 2>&1`
	echo "$time"
	sum=$(echo $sum + $time | bc -l);
done

avg=$(echo $sum / $loop_variable | bc -l);
avg=`printf "%.1f" $avg`
echo "Average time: $avg"

# 5 4 3

echo "--- 5 4 3 ---"

sum=0
for (( i=0; i<loop_variable; i++ ))
do
	time=`./ballAlg 5 4 3 2>&1`
	echo "$time"
	sum=$(echo $sum + $time | bc -l);
done

avg=$(echo $sum / $loop_variable | bc -l);
avg=`printf "%.1f" $avg`
echo "Average time: $avg"

# 2 6 8

echo "--- 2 6 8 ---"

sum=0
for (( i=0; i<loop_variable; i++ ))
do
	time=`./ballAlg 2 6 8 2>&1`
	echo "$time"
	sum=$(echo $sum + $time | bc -l);
done

avg=$(echo $sum / $loop_variable | bc -l);
avg=`printf "%.1f" $avg`
echo "Average time: $avg"