#!/bin/bash
#SBATCH --time=24:00:00                # Ensure the time is sufficient for each job
#SBATCH --mem=16G                      # Ensure the memory is sufficient for each job
#SBATCH --nodes=1                      # Number of nodes per job
#SBATCH --cpus-per-task=2              # Number of CPUs per task

# Load the necessary module
module load ROOT

# Get the pthardbin value from the command-line arguments
pthardbin=$1

# Navigate to the directory containing your ROOT script
cd /home/ar2545/RooUnfold/examples/


# Run your ROOT script with the pthardbin value
root -l -b -q 'Unfolding3dAnalysis.cxx.cxx("",'$pthardbin')'

# You might need to update the above line to match how your script handles the pthardbin value.
