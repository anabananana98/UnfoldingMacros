# UnfoldingMacros
Unfolding scripts 

To run the 3D unfolding script simply go into /home/ar2545/RooUnfold on GRACE and run 

```bash
sbatch trythreeDunf.sh
```

This script runs over each pTHardBin separately. It calls the macro named Unfolding3dAnalysis.cxx and uses the source script called run3Dunfolding.sh
These files only build the response matrix because eventually to unfold we need the full response.

For performing the actual unfolding I have a script called “PerformUnfolding.cxx”. The script /home/ar2545/RooUnfold/doUnfold.sh can be sourced to do this and the unfolding jobs can be launched by running 
```bash
source launchdoUnfold.sh
```

CompareResults.C and CompareResultsNIter.C compares unfolded results vs binbybin corrected results.
