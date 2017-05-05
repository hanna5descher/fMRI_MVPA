clear;
curDir = pwd;

%% SETTINGS

nameAnalysis = 'MVPA_ObjectiveCueValueRegression_WithinSubject_NativeSpace_perPhaseGLM';
nameRegressor = 'MVPA_CueValue_pmodFeedbackRT_perRun_NativeSpace_perPhase'; % this should match the folder name of GLM analysis

% define searchlight file to use (size of the sphere)
nameIdx = 'searchIdx4'; % size of searchlight sphere, e.g., r = 2 voxels for searchIdx2.dat
fileID = 'beta'; % beta or Tmap
dirIdxNative = fullfile(curDir,'searchIdx','NativeSpace'); % where searchlight location files are stored
infoIdx = spm_vol(fullfile(dirIdxNative,'fSAT01_spmT_0001.hdr'));
fwhm = 5; % for smoothing

% set directories
dirData = fullfile(curDir,'dataFMRI');
dirFeatures = fullfile(curDir,'MVPA','Features');
dirRegInfo = fullfile(dirFeatures,nameRegressor);
if ~exist(dirRegInfo,'dir')
    error('warning: there is no matching folder for this analysis');
else
    load(fullfile(dirRegInfo,['reg' nameRegressor])); % load regressor and contrast information
end
dirSaveFeatures = fullfile(dirRegInfo,[fileID '_' nameIdx]); % where features will be saved
dirAccuracy = fullfile(curDir,'MVPA','AccuracyMap');
dirSaveAccuracy = fullfile(dirAccuracy,[nameAnalysis '_smoothed_s' num2str(fwhm)],[fileID '_' nameIdx]);% where accuracy maps will be saved
if ~exist(dirSaveAccuracy,'dir')
    mkdir(dirSaveAccuracy);
end

%% LOAD FEATURES and RUN MVPA

nameAccuracyMap = {'NP_CueRegression','TP_CueRegression'};

for iS = 1 : scanP.nSubj
    dirIdx{iS} = fullfile(dirIdxNative,[scanP.subjID{iS} '_' nameIdx '.dat']); % full path to the searchlight file
    dirHdr{iS} = fullfile(dirIdxNative,['f' scanP.subjID{iS} '_spmT_0001.hdr']);
    
    tempF = load(fullfile(dirSaveFeatures,scanP.subjID{iS}),'features'); % load feature matrix
    
    for iP = 1 : scanP.nPhase
        F = tempF.features{iP};
        label = [0.8; 0.6; 0.4; 0.2]; % label for each feature
        
        disp(['===== Running searchlight: ' [scanP.subjID{iS} '_' nameAccuracyMap{iP}] ' =====']);
        mvpaRunSearchlightRegressionWithinSubjects(dirIdx{iS},infoIdx.dim,F,label,scanP.nRunPhase,dirHdr{iS},dirSaveAccuracy,[nameAccuracyMap{iP} '_' scanP.subjID{iS}]);
       
        clear F label
    end
    
    clear tempF 
end

save(fullfile(dirSaveAccuracy,['info' nameAnalysis]),'nameRegressor','dirHdr','dirIdx','dirRegInfo','dirSaveFeatures','dirSaveAccuracy','nameAccuracyMap');

%cd(dirSaveAccuracy);

%% NORMALIZE ACCURACY & MSE MAPS INTO THE MNI TEMPLATE
% load(fullfile(dirSaveAccuracy,['info' nameAnalysis]));
 
spm('defaults', 'FMRI');
for iS = 1 : scanP.nSubj
    % NORMALIZE ACCURACYINTO THE MNI TEMPLATE
    dirAnat = fullfile(dirData,scanP.subjID{iS},scanP.biacID{iS},'Anat'); % T1 anatomical
    
    % locate deformation field (created as y*.nii)
    nameDeformField = {fullfile(dirAnat,['y_' scanP.subjID{iS} '_anat0001.nii'])};
    
    % locate images
    for iP = 1 : scanP.nPhase
        nameAccuracy(iP,:) = {[fullfile(dirSaveAccuracy,['cv' nameAccuracyMap{iP} '_' scanP.subjID{iS} '.img']) ',1']};
    end   
    matlabbatch{1}.spm.spatial.normalise.write.subj.def = nameDeformField;
    matlabbatch{1}.spm.spatial.normalise.write.subj.resample = nameAccuracy;
    matlabbatch{1}.spm.spatial.normalise.write.woptions.bb = [-78 -112 -70; 78 76 85];
    matlabbatch{1}.spm.spatial.normalise.write.woptions.vox = [3,3,3]; % resample as the original resolution
    matlabbatch{1}.spm.spatial.normalise.write.woptions.interp = 4;
    matlabbatch{1}.spm.spatial.normalise.write.woptions.prefix = 'w';
    
    % SMOOTHING
    % locate normalized functional images
    for iP = 1 : scanP.nPhase
        nameNormAccuracy(iP,:) = {[fullfile(dirSaveAccuracy,['wcv' nameAccuracyMap{iP} '_' scanP.subjID{iS} '.img']) ',1']};
    end    
    matlabbatch{2}.spm.spatial.smooth.data = nameNormAccuracy;
    matlabbatch{2}.spm.spatial.smooth.fwhm = repmat(fwhm,1,3); % smoothing kernel size = 6
    matlabbatch{2}.spm.spatial.smooth.dtype = 0;
    matlabbatch{2}.spm.spatial.smooth.im = 0;
    matlabbatch{2}.spm.spatial.smooth.prefix = 's';
    
    spm_jobman('serial',matlabbatch);
    
    clear matlabbatch nameDeformField nameAccuracy 
end

%% 2ND LEVEL PAIRED T-TEST
% load(fullfile(dirSaveAccuracy,['info' nameAnalysis]));
dirSaveGroupAnalysis = fullfile(dirSaveAccuracy,'2ndLevel_PairedTTest'); % where design matrix for the 2nd level analysis will be saved
if ~exist(dirSaveGroupAnalysis,'dir')
    mkdir(dirSaveGroupAnalysis);
end

matlabbatch{1}.spm.stats.factorial_design.dir = {dirSaveGroupAnalysis};
for iS = 1 : scanP.nSubj
    for iP = 1 : scanP.nPhase
        matlabbatch{1}.spm.stats.factorial_design.des.pt.pair(iS).scans(iP,:) = {fullfile(dirSaveAccuracy,['swcv' nameAccuracyMap{iP} '_' scanP.subjID{iS} '.img,1'])};
    end    
end

matlabbatch{1}.spm.stats.factorial_design.des.pt.gmsca = 0; % grand mean scaling: not applicable for fmri data
matlabbatch{1}.spm.stats.factorial_design.des.pt.ancova = 0; % ancova: not applicable for fmri data
matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {}); % covariates
matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {}); % multiple covariates
matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1; % threhsold masking: none
matlabbatch{1}.spm.stats.factorial_design.masking.im = 1; % implicit masking: yes (NaNs treated as missing)
matlabbatch{1}.spm.stats.factorial_design.masking.em = {''}; % explicit masking: no
matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1; % global calculation: omit
matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1; % grand mean scaling: no
matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1; % global normalization: no

% estimate
matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep('Factorial design specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;

% contrast
matlabbatch{3}.spm.stats.con.spmmat(1) = cfg_dep('Model estimation: SPM.mat File', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{3}.spm.stats.con.consess{1}.tcon.name = 'PairedT';
matlabbatch{3}.spm.stats.con.consess{1}.tcon.weights = [1 -1];
matlabbatch{3}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
matlabbatch{3}.spm.stats.con.delete = 0;

spm('defaults', 'FMRI');
spm_jobman('serial',matlabbatch)

%% 2ND LEVEL ONE-SAMPLE T-TEST
%load(fullfile(dirSaveAccuracy,['info' nameAnalysis]));
clear matlabbatch
nameGroupFolder = {'2ndLevel_NP','2ndLevel_TP'};
prefixGroup = {'swcvNP', 'swcvTP'};

spm('defaults', 'FMRI');
for iG = 1 : length(nameGroupFolder)
    dirSaveGroupAnalysis = fullfile(dirSaveAccuracy,nameGroupFolder{iG}); % where design matrix for the 2nd level analysis will be saved
    if ~exist(dirSaveGroupAnalysis,'dir')
        mkdir(dirSaveGroupAnalysis);
    end
    
    matlabbatch{1}.spm.stats.factorial_design.dir = {dirSaveGroupAnalysis};
    for iS = 1 : scanP.nSubj
        matlabbatch{1}.spm.stats.factorial_design.des.t1.scans(iS,:) = {fullfile(dirSaveAccuracy,[prefixGroup{iG} '_CueRegression_' scanP.subjID{iS} '.img,1'])};
    end
    
    matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
    matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
    matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
    matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
    matlabbatch{1}.spm.stats.factorial_design.masking.em = {''};
    matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
    matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
    matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;
    matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep('Factorial design specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
    matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
    matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;
    matlabbatch{3}.spm.stats.con.spmmat(1) = cfg_dep('Model estimation: SPM.mat File', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
    matlabbatch{3}.spm.stats.con.consess{1}.tcon.name = 'T';
    matlabbatch{3}.spm.stats.con.consess{1}.tcon.weights = 1;
    matlabbatch{3}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
    matlabbatch{3}.spm.stats.con.delete = 0;
    
    spm_jobman('serial',matlabbatch)
    
    clear dirSaveGroupAnalysis matlabbatch
end

cd(dirSaveAccuracy)


