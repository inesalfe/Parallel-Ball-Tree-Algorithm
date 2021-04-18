#!/bin/sh

loop_variable=$5
os="1"

# 2 5 1

echo "--- 2 5 0 ---"

./ballAlg 2 5 "$os" 2>&1
output=`./ballQuery pts.txt 3 1 2>&1`

if [[ "$output" =~ "2.777747 5.539700" ]]; then
  echo "Correct output"
fi

# 2 8 1

echo "--- 2 8 0 ---"

./ballAlg 2 8 "$os" 2>&1
output=`./ballQuery pts.txt 8 8 2>&1`

if [[ "$output" =~ "7.830992 7.984400" ]]; then
  echo "Correct output"
fi

# 20 1000000 0

echo "--- 20 1000000 0 ---"

./ballAlg 20 1000000 "$os" 2>&1
output=`./ballQuery pts.txt 1 2 3 4 5 6 7 8 9 1 2 3 4 5 6 7 8 9 1 2 2>&1`

if [[ "$output" =~ "0.809638 0.609105 2.107746 2.097437 5.055056 7.446720 5.392867 8.758850 9.246286 3.082144 0.496353 2.529218 0.473570 3.954821 6.931734 5.535411 8.118067 9.333055 0.469908 2.162746" ]]; then
  echo "Correct output"
fi

# 3 5000000 0

echo "--- 3 5000000 0 ---"

./ballAlg 3 5000000 "$os" 2>&1
output=`./ballQuery pts.txt 4 5 6 2>&1`

if [[ "$output" =~ "3.979046 5.032039 6.011886" ]]; then
  echo "Correct output"
fi

# 4 10000000 0

echo "--- 4 10000000 0 ---"

./ballAlg 4 10000000 "$os" 2>&1
output=`./ballQuery pts.txt 2 4 6 8 2>&1`

if [[ "$output" =~ "1.996719 4.012344 5.988101 8.081113" ]]; then
  echo "Correct output"
fi

# 3 20000000 0

echo "--- 3 20000000 0 ---"

./ballAlg 3 20000000 "$os" 2>&1
output=`./ballQuery pts.txt 1 5 9 2>&1`

if [[ "$output" =~ "1.003042 4.986528 9.010856" ]]; then
  echo "Correct output"
fi

# 4 20000000 0

echo "--- 4 20000000 0 ---"

./ballAlg 4 20000000 "$os" 2>&1
output=`./ballQuery pts.txt 8 6 4 2 2>&1`

if [[ "$output" =~ "7.939939 5.934679 3.951869 1.930474" ]]; then
  echo "Correct output"
fi

# 3 3 3

echo "--- 3 3 3 ---"

./ballAlg 3 3 3 2>&1
output=`./ballQuery pts.txt 5 5 5 2>&1`

if [[ "$output" =~ "5.613802 2.249833 3.930918" ]]; then
  echo "Correct output"
fi

# 5 4 3

echo "--- 5 4 3 ---"

./ballAlg 5 4 3 2>&1
output=`./ballQuery pts.txt 9 7 5 3 1 2>&1`

if [[ "$output" =~ "5.613802 2.249833 3.930918 4.439384 2.850413" ]]; then
  echo "Correct output"
fi

# 2 6 8

echo "--- 2 6 8 ---"

./ballAlg 2 6 8 2>&1
output=`./ballQuery pts.txt 3 7 2>&1`

if [[ "$output" =~ "3.527607 7.895896" ]]; then
  echo "Correct output"
fi

