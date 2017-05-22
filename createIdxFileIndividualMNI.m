clear;
%addpath('\dataBehavioral');
%load('dataExtractBehavior3','scanP'); % load scan parameters
dirData = fullfile(pwd,'dataFMRI');
curDir = pwd;
dirSave = fullfile(curDir,'IndividualMNI');

nameAnalysis = 'MVPA_CueValue_pmodFeedbackRT_perRun_MNIs5wa_perPhase';

for iS = 1 : scanP.nSubj
    dirAnat = fullfile(dirData,scanP.subjID{iS},scanP.biacID{iS},'Anat'); % T1 anatomical
    dirFunc = fullfile(dirData,scanP.subjID{iS},scanP.biacID{iS},'Func'); % functional scans
     
    % normalize grey matter mask
    nameDeformField = {fullfile(dirAnat,['y_' scanP.subjID{iS} '_anat0001.nii'])};
    nameGM = {[fullfile(dirAnat,['c1' scanP.subjID{iS} '_anat0001.nii']) ',1']}; % maps full directory of T1 image
   
    matlabbatch{1}.spm.spatial.normalise.write.subj.def = nameDeformField;
    matlabbatch{1}.spm.spatial.normalise.write.subj.resample = nameGM;
    matlabbatch{1}.spm.spatial.normalise.write.woptions.bb = [-78 -112 -70; 78 76 85];
    matlabbatch{1}.spm.spatial.normalise.write.woptions.vox = [1 1 1]; % resample as the original resolution
    matlabbatch{1}.spm.spatial.normalise.write.woptions.interp = 4;
    matlabbatch{1}.spm.spatial.normalise.write.woptions.prefix = 'w';
    
    spm('defaults', 'FMRI');
    spm_jobman('serial',matlabbatch);
    clear matlabbatch
    
    % resample the mask to match functional space
    nameGM_MNI = fullfile(dirAnat,['wc1' scanP.subjID{iS} '_anat0001.nii']);
    nameGM_MNI = fullfile(curDir,'c1T1_Original.nii'); % to segment template GM mask
    nameTargetSpace = fullfile(dirFunc,'Run1',['wa' scanP.subjID{iS} '_Run1_0005.img']);
    nameImages = {nameTargetSpace, nameGM_MNI};
    optReslice.which = [1 0]; % images 2...n (don't reslice the sample image)
    optReslice.interp = 0; % nearest neighbor
    optReslice.wrap = [0 0 0];
    optReslice.mask = 1;
    optReslice.prefix = 'r';
    spm_reslice(nameImages,optReslice);
    
    % create searchlight index file    
    % if needed, convert .nii file to .hdr/.img files; GenerateIdxFile processes .img file only.
    nameGMResampled = fullfile(dirAnat,['frwc1' scanP.subjID{iS} '_anat0001.img']);
    %nameGMResampled = fullfile(curDir,'frc1T1_Original.img');
    infoHdr = spm_vol(nameGMResampled);
    threshold = 1;
    r = 3; % searchlight radus (e.g., r = 3 voxels)
    saveName = fullfile(dirSave,[scanP.subjID{iS} '_searchIdx' num2str(r) '.dat']);
    GenerateIdxFile(nameGMResampled, infoHdr.dim, saveName, threshold, r);

    % copy spmT files to create header files
    dirDataT = fullfile(dirData,scanP.subjID{iS},'Analysis',nameAnalysis,'Phase1_NP','beta_0001.nii'); % location of the header file
    copyfile(dirDataT,fullfile(dirSave,[scanP.subjID{iS} '_spm_beta_0001.nii']));
end
