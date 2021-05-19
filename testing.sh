#!/bin/sh

os="0"
proc="8"

# 20 1000000 0

./ballQuery pts1.txt 1 2 3 4 5 6 7 8 9 1 2 3 4 5 6 7 8 9 1 2

echo "--- 20 1000000 0 ---"

output=`./ballQuery <(srun --ntasks-per-node=1 -n "$proc" ballAlg 20 1000000 "$os") 1 2 3 4 5 6 7 8 9 1 2 3 4 5 6 7 8 9 1 2 2>&1`

if [[ "$output" =~ "0.809638 0.609105 2.107746 2.097437 5.055056 7.446720 5.392867 8.758850 9.246286 3.082144 0.496353 2.529218 0.473570 3.954821 6.931734 5.535411 8.118067 9.333055 0.469908 2.162746" ]]; then
  echo "Correct output"
else
  echo "Incorrect output"
fi

# 3 5000000 0

./ballQuery pts2.txt 4 5 6

echo "--- 3 5000000 0 ---"

output=`./ballQuery <(srun --ntasks-per-node=1 -n "$proc" ballAlg 3 5000000 "$os") 4 5 6 2>&1`

if [[ "$output" =~ "3.979046 5.032039 6.011886" ]]; then
  echo "Correct output"
else 
  echo "Incorrect output"
fi

# 4 10000000 0

./ballQuery pts3.txt 2 4 6 8

echo "--- 4 10000000 0 ---"

output=`./ballQuery <(srun --ntasks-per-node=1 -n "$proc" ballAlg 4 10000000 "$os") 2 4 6 8 2>&1`

if [[ "$output" =~ "1.996719 4.012344 5.988101 8.081113" ]]; then
  echo "Correct output"
else 
  echo "Incorrect output"
fi

# 3 20000000 0

./ballQuery pts4.txt 1 5 9

echo "--- 3 20000000 0 ---"

output=`./ballQuery <(srun --ntasks-per-node=1 -n "$proc" ballAlg 3 20000000 "$os") 1 5 9 2>&1`

if [[ "$output" =~ "1.003042 4.986528 9.010856" ]]; then
  echo "Correct output"
else 
  echo "Incorrect output"
fi

# 4 20000000 0

./ballQuery pts5.txt 8 6 4 2

echo "--- 4 20000000 0 ---"

output=`./ballQuery <(srun --ntasks-per-node=1 -n "$proc" ballAlg 4 20000000 "$os") 8 6 4 2 2>&1`

if [[ "$output" =~ "7.939939 5.934679 3.951869 1.930474" ]]; then
  echo "Correct output"
else 
  echo "Incorrect output"
fi
