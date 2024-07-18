#!/bin/bash
#SBATCH --job-name=ThreeDUnfolding
#SBATCH --time=24:00:00
#SBATCH --mem=16G
#SBATCH --nodes=2
#SBATCH --cpus-per-task=2

# Load the necessary module
module load ROOT

# Define the array of pthardbin values
pthardbins=(2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20)

# Iterate over each pthardbin value and submit it as a separate job
for pthardbin in "${pthardbins[@]}"
do
    sbatch --job-name="Unfolding3dAnalysis_${pthardbin}" \
           --output="job_output_%j_${pthardbin}.out" \
           --error="job_error_%j_${pthardbin}.err" \
           /home/ar2545/RooUnfold/run3Dunfolding.sh "$pthardbin"
done
