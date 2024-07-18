#!/bin/bash
#SBATCH --job-name=“ThreeDUnfolding”
#SBATCH --time=24:00:00
#SBATCH --mem=16G
#SBATCH --nodes=2
#SBATCH --cpus-per-task=2

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






#root -l -b -q '/home/ar2545/RooUnfold/examples/Unfolding3dAnalysis.cxx("", 2)'

#root -l -b -q '/home/ar2545/RooUnfold/examples/Unfolding3dAnalysis.cxx("", 3)'
#root -l -b -q '/home/ar2545/RooUnfold/examples/Unfolding3dAnalysis.cxx("", 4)'
#root -l -b -q '/home/ar2545/RooUnfold/examples/Unfolding3dAnalysis.cxx("", 5)'
#root -l -b -q '/home/ar2545/RooUnfold/examples/Unfolding3dAnalysis.cxx("", 6)'
#root -l -b -q '/home/ar2545/RooUnfold/examples/Unfolding3dAnalysis.cxx("", 7)'


# Define the tag argument for your ROOT macro
#tag=""  # Replace with your tag argument

# Array of pthardbin values to iterate over
#pthardbins=(2 3)  # Example list of pthardbin values

# Iterate over each pthardbin value and run the Bash script
#for pthardbin in "${pthardbins[@]}"
#do
#   /home/ar2545/RooUnfold/examples/Unfolding3dAnalysis.cxx() tag $pthardbin &
#done



