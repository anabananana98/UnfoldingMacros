#!/bin/bash
#SBATCH --job-name=Unf3D
#SBATCH --output=Unf3D_output_%j.txt
#SBATCH --error=Unf3D_error_%j.txt
#SBATCH --ntasks=1
#SBATCH --time=24:00:00
#SBATCH --mem=16G

# Load necessary modules
module load ROOT

# Assign command-line arguments to variables
RESPONSE_FILE=$1
DATA_FILE=$2
N_ITER=$3

# Check if the correct number of arguments are provided
if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <response_file> <data_file> <n_iterations>"
    exit 1
fi

# Navigate to the directory containing your ROOT script
cd /home/ar2545/RooUnfold/examples/

# Run the ROOT macro
root -l -b -q "Perform3DUnfolding.cxx(\"${RESPONSE_FILE}\", \"${DATA_FILE}\", ${N_ITER})"
