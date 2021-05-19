#!/bin/sh
#SBATCH --job-name=ballAlg_i_r
#SBATCH --error=stderr
#SBATCH --ntasks=16
#SBATCH --ntasks-per-node=1

### 12 563476 3

echo "--- 12 563476 3 ---"

output=`./ballQuery <(srun ballAlg 12 563476 3 &) 7 3 4 5 6 1 3 6 2 5 1 2 2>&1`

if [[ "$output" =~ "5.892031 3.743222 5.078488 3.908186 5.590736 2.332848 2.638832 5.872337 2.048928 6.505624 1.709915 1.254254" ]]; then
  echo "Correct output"
else
  echo "Incorrect output"
fi

# ### 7 750233 2

# echo "--- 7 750233 2 ---"

# output=`./ballQuery <(srun 7 750233 2) 4 5 6 7 1 3 2 2>&1`

# if [[ "$output" =~ "4.245852 5.516109 5.973159 7.386220 0.625611 2.495845 1.520084" ]]; then
#   echo "Correct output"
# else 
#   echo "Incorrect output"
# fi

# ### 4 600000 0

# echo "--- 4 600000 0 ---"

# output=`./ballQuery <(srun ballAlg 4 600000 "$os") 2 4 6 8 2>&1`

# if [[ "$output" =~ "1.874237 3.944711 5.964251 7.886127" ]]; then
#   echo "Correct output"
# else 
#   echo "Incorrect output"
# fi

# ### 20 432167 4

# echo "--- 20 432167 4 ---"

# output=`./ballQuery <(srun ballAlg 20 432167 4) 1 5 9 5 4 7 3 4 5 6 7 8 2 3 4 1 3 4 2 5 2>&1`

# if [[ "$output" =~ "3.292571 3.887630 7.860798 4.393064 3.032862 3.084203 0.604753 5.289977 3.913724 7.467245 6.843754 6.489080 3.223100 0.937685 4.450096 1.946901 1.472640 4.699993 2.725623 3.598386" ]]; then
#   echo "Correct output"
# else 
#   echo "Incorrect output"
# fi

# ### 3 250000 5

# echo "--- 3 250000 5 ---"

# output=`./ballQuery <(srun ballAlg 3 250000 5) 8 6 4 2>&1`

# if [[ "$output" =~ "7.979664 5.983951 3.969156" ]]; then
#   echo "Correct output"
# else 
#   echo "Incorrect output"
# fi
