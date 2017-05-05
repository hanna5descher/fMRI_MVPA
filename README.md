# fMRI_MVPA
Whole-Brain Searchlight Multi-Voxel Pattern Analysis (MVPA) using Support Vector Machine Classification &amp; Regression

Included MATLAB codes use SPM12 (or SPM8; http://www.fil.ion.ucl.ac.uk/spm/) for functional magnetic resonance imaging (fMRI) data preprocessing and the general linear model (GLM) analysis.

## Analysis Pipeline
1. fMRI data preprocessing: PreprocessingBatchSPM12.m (or PreprocessingSPM8.m for older SPM versions)
   
   Functional image realignment (and reslice, if running MVPA in a subject-specific native space), slice timing correction, coregistration, T1 segmentation, normalizatio (to MNI space), and smoothing (optional)

2. Run GLM
   
   To construct inputs (or features) for MVPA, perform GLM analysis to extract beta or T maps per condition
   * Native space: Use resliced functional images (I prefer to use unsmoothed data)
   * MNI space: Use normalized functional images (again, using unsmoothed data is more conventional)

3. Build the whole-brain searchlight using a gray matter mask: createIdxFileIndividualMNI.m or createIdxFileNativeSpace.m
   
   To perform searchlight MVPA (Kriegeskorte et al. 2006; Haynes et al. 2007), generate a whole-brain searchlight file using a gray matter mask from the T1 segmentation. The searchlight is a spherical cluster with a radius of 3-5 voxels.
   * Resample gray matter masks (either from subject's T1 or from the SPM T1 template) to match the functional image space
   * Create a searchlight index file (.dat): GenerateIdxFile.m
   * Note. The function "GenerateIdxFile" uses a brute force method (lots of for loops) to search through the entire brain...(could be enhanced later, but it works!)

4. Feature generation
   
   Using beta (mvpaMakeFeaturesBeta.m) or T (mvpaMakeFeaturesT.m) maps (#2) and the searchlight file (#3), extract features for MVPA.

5. Run SVM Classification (mvpaRunSearchlightWithinSubjects.m) or Regression (mvpaRunSearchlightRegressionWithinSubjects.m)

6. Generate a group-level average accuracy map to detect statistically significant clusters (= searchlights)

* runMVPA_SVM_Regression.m is a sample script that covers the procedure decribed in 3-6 for SVM regression. SVM classification uses the same procedure by switching out the "label" variable to binary (0 and 1) and using the function, "mvpaRunSearchlightWithinSubjects".
