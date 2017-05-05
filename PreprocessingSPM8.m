clear; close all;
% Preprocessing Script for SPM8 (or prior)

load('dataAnalyzeBehavior','subjID','biacID');
for iS = 1 : length(subjID)
scanInfo.subjID = subjID{iS}; %'SAT32';
scanInfo.biacID = biacID{iS}; %'20150515_19775'; % BIAC ID for the functional scan
scanInfo.dirData = fullfile(pwd,'Data'); % fMRI data directory
scanInfo.dirSaveData = fullfile(scanInfo.dirData,scanInfo.subjID,scanInfo.biacID);
scanInfo.dirAnat = fullfile(scanInfo.dirData,scanInfo.subjID,scanInfo.biacID,'Anat'); % location of T1 anatomical
scanInfo.dirFunc = fullfile(scanInfo.dirData,scanInfo.subjID,scanInfo.biacID,'Func'); % location of functional scans
scanInfo.TR = 2; % TR in sec.
scanInfo.nTRDisdaq = 4; % number of TRs to discard
scanInfo.nTRRun = [165 165 165 165 135 135 135 135 215 150]; % number of TRs for each session (not including disdaq)
scanInfo.nRun = length(scanInfo.nTRRun); % total number of runs
scanInfo.nSlice = 42; % number of slices (for slice timing correction)
scanInfo.resolution = 3; % spatial resolution of normalized fmri data: should use the same resolution as the functional image (prefix w3 = 3; w = 2 <- this is no bueno)
scanInfo.smoothingKernelSize = 8; % size of Gaussian kernel for smoothing (e.g., 8)
scanInfo.MVPAres = 3; % resampling for MVPA; use the actual voxel size (e.g., 3)
cd(scanInfo.dirSaveData);

% read in raw functional files
scanInfo.nameFunc = cell(1,scanInfo.nRun);
for nR = 1 : scanInfo.nRun
    scanInfo.dirFuncRuns{nR} = fullfile(scanInfo.dirFunc,['Run' num2str(nR)]); % maps location of files for each run
    nameFiles = spm_select('List',scanInfo.dirFuncRuns{nR},'.img'); % detects .img file names
    for nT = 1 : scanInfo.nTRRun(nR)
        scanInfo.nameFunc{1,nR}(nT,:) = fullfile(scanInfo.dirFuncRuns{nR},nameFiles(nT+scanInfo.nTRDisdaq,:)); % maps full directory name for each funcitonal file and discards disdaq
    end
end

% read in T1 anantomical file
nameFile = spm_select('List',scanInfo.dirAnat,'.img');
scanInfo.nameAnat = fullfile(scanInfo.dirAnat,nameFile(1,:)); % maps full directory of T1 image

spm('defaults', 'FMRI')
spm_figure('Create','Graphics');

disp(['=========== Begin Preprocessing for ' scanInfo.subjID ' ==========']);

%% REALIGNMENT (ESTIMATE & RESLICE)
% Motion correction
% Do realignment and then slice timing 
% see http://imaging.mrc-cbu.cam.ac.uk/imaging/SliceTimingRikEmail for a discussion on this issue.
% Realignment creates text files with 6 columns: [translations in mm (”right”,“forward”,“up”), estimated rotations in rad (”pitch”,“roll”,“yaw”)]
% Performing this step will also create a mean image from your set of functional images, 
% which can later be used to match functional to structural images.
disp('Step 1: Realignment (Estimate & Reslice)');

% Realignment
paramRealignEst.quality = 0.9;
paramRealignEst.sep = 4;
paramRealignEst.fwhm = 5;
paramRealignEst.rtm = 1; % register to mean (don't change it!)
paramRealignEst.interp = 2;
paramRealignEst.wrap = [0 0 0];
paramRealignEst.graphics = 1; % outputs headmotion graph
spm_realign(scanInfo.nameFunc,paramRealignEst);

% Reslice
paramRealignWrite.which = [0 1]; % mean images only
paramRealignWrite.interp = 4;
paramRealignWrite.wrap = [0 0 0];
paramRealignWrite.mask = 0; % no masking
paramRealignWrite.prefix = 'r';
spm_reslice(scanInfo.nameFunc,paramRealignWrite);

disp('=========== Step 1 Done ==========');

% SLICE TIMING
disp('Step 2: Slice Timing');

scanInfo.sliceOrder = [1:2:scanInfo.nSlice, 2:2:scanInfo.nSlice]; % slice order for Slice Timing
scanInfo.refSlice = 1;
scanInfo.TA = scanInfo.TR - (scanInfo.TR/scanInfo.nSlice);
scanInfo.timing(1) = scanInfo.TA / (scanInfo.nSlice-1);
scanInfo.timing(2) = scanInfo.TR - scanInfo.TA;
scanInfo.prefixSliceTiming = 'a';
spm_slice_timing(scanInfo.nameFunc,scanInfo.sliceOrder,scanInfo.refSlice,scanInfo.timing,scanInfo.prefixSliceTiming);

disp('=========== Step 2 Done ==========');

% COREGISTER (ESTIMATE)
% Overlay structural and functional images: Link functional scans to anatomical scan
% REFERENCE IMAGE: read in mean functional image (remains stationary)
% SOURCE IMAGE: read in anatomical T1 image (image to be jiggled)
% disp('Step 3: Coregister (Estimate)');

locateMean = dir(fullfile(scanInfo.dirFunc,'Run1','mean*.img')); % mean file is created within the Run1 folder
scanInfo.nameFuncMean = fullfile(scanInfo.dirFunc,'Run1',locateMean.name);
paramCoreg.cost_fun = 'nmi'; % Normalised Mutual Information
paramCoreg.sep = [4 2];
paramCoreg.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
paramCoreg.fwhm = [7 7];
paramCoreg.graphics = 1;

x  = spm_coreg(spm_vol(scanInfo.nameFuncMean),spm_vol(scanInfo.nameAnat),paramCoreg);
M  = inv(spm_matrix(x));

% coregister anatomical to mean functional
MM = spm_get_space(deblank(scanInfo.nameAnat));
spm_get_space(deblank(scanInfo.nameAnat), M*MM);

disp('=========== Step 3 Done ==========');


%% NORMALIZE (ESTIMATE & REWRITE)
% Warp images to fit to a standard template brain (MNI)
disp('Step 4: Normalize (Estimate & Rewrite)');

% Normalize anatomical to MNI template
templateName = 'C:\Users\ho20\Documents\MATLAB\spm8\templates\T1.nii,1';
paramNormEst.smosrc = 8;
paramNormEst.smoref = 0;
paramNormEst.regtype = 'mni';
paramNormEst.cutoff = 25;
paramNormEst.nits = 16;
paramNormEst.reg = 1; %regulariazation factor. Default is 1. Increase the value exponentially if your normalization seems distorted.
sn = spm_normalise(templateName,scanInfo.nameAnat,fullfile(scanInfo.dirAnat, [scanInfo.subjID '_anat0001_sn.mat']), '', '',paramNormEst); 

% Rewrite EPI images (use corrected images with prefix 'a')
disp('Writing normalized EPIs...');
scanInfo.nameFuncA = cell(1,scanInfo.nRun);
prefix = sprintf('^%s.*\\.img', 'a'); % input slice time corrected images
for nR = 1 : scanInfo.nRun
    nameFiles = spm_select('List',scanInfo.dirFuncRuns{nR},prefix); % detects a*.img file names
    for nT = 1 : scanInfo.nTRRun(nR)
        nameFuncA{1,nR}(nT,:) = fullfile(scanInfo.dirFuncRuns{nR},nameFiles(nT,:)); % maps full directory name for each funcitonal file and discards disdaq
    end
end
scanInfo.nameFuncA = strvcat(nameFuncA); % gets rid of cells and arranges names into a single matrix
paramNormWrite.preserve = 0;
paramNormWrite.bb = [-78 -112 -50; 78   76   85];
paramNormWrite.vox = repmat(scanInfo.resolution,1,3);
paramNormWrite.interp = 1;
paramNormWrite.wrap = [0 0 0];
paramNormWrite.prefix = 'w';
scanInfo.prefixNormalize = paramNormWrite.prefix;
for i = 1 : size(scanInfo.nameFuncA,1)
    spm_write_sn(scanInfo.nameFuncA(i,:),sn,paramNormWrite);
end

% T1 Anatomical Rewrite 
% This is to view the quality of normalization using T1 (better resolution than functional scans)
disp('Writing T1...');
paramNormWriteT1 = paramNormWrite;
paramNormWriteT1.vox = [1 1 1];
spm_write_sn(scanInfo.nameAnat,sn,paramNormWriteT1);

% % PreProcessNormalizeResample (this is for MVPA, rewrite with an actual EPI voxel size)
% disp('Resampling...');
% paramNormResample.preserve = 0;
% paramNormResample.bb = [-78 -112 -50; 78   76   85];
% paramNormResample.vox = repmat(scanInfo.MVPAres,1,3);
% paramNormResample.interp = 1;
% paramNormResample.wrap = [0 0 0];
% paramNormResample.prefix = 'rw';
% scanInfo.prefixNormResample = paramNormResample.prefix;
% for i = 1 : size(scanInfo.nameFuncA,1)
%     spm_write_sn(scanInfo.nameFuncA(i,:),sn,paramNormResample);
% end

disp('=========== Step 4 Done ==========');

%% SMOOTHING
% To increase signal-to-noise ratio
% use normalized images with prefix 'wa'
disp('Step 5: Smoothing');
scanInfo.prefixSmooth = 's';
nameFuncWA = cell(1,scanInfo.nRun);
nameFuncOut = cell(1,scanInfo.nRun);
prefix = sprintf('^%s.*\\.img', 'wa'); % input normalized images
for nR = 1 : scanInfo.nRun
    nameFiles = spm_select('List',scanInfo.dirFuncRuns{nR},prefix); % detects a*.img file names
    for nT = 1 : scanInfo.nTRRun(nR)
        nameFuncWA{1,nR}(nT,:) = fullfile(scanInfo.dirFuncRuns{nR},nameFiles(nT,:)); % maps full directory name for each funcitonal file and discards disdaq
        nameFuncOut{1,nR}(nT,:) = fullfile(scanInfo.dirFuncRuns{nR},[scanInfo.prefixSmooth nameFiles(nT,:)]);
    end
end
scanInfo.nameFuncWA = strvcat(nameFuncWA);
nameFilesSmooth = strvcat(nameFuncOut); % name of smoothed output files
paramSmooth.fwhm = repmat(scanInfo.smoothingKernelSize,1,3);
paramSmooth.dtype = 0;
for i = 1 : size(scanInfo.nameFuncWA,1)
    spm_smooth(scanInfo.nameFuncWA(i,:),nameFilesSmooth(i,:),paramSmooth.fwhm,paramSmooth.dtype);
end

disp('=========== Step 5 Done ==========');

%% SAVE PREPROCESSING INFORMATION
%save([scanInfo.subjID '_dataPreprocess.mat'],'scanInfo','paramRealignEst','paramRealignWrite','paramCoreg','paramNormEst','paramNormWrite','paramNormResample','paramSmooth');
save([scanInfo.subjID '_dataPreprocessRedoVox3.mat'],'scanInfo','paramNormEst','paramNormWrite','paramNormWriteT1','paramSmooth');

disp('=========== Preprocessing Completed! ==========');

clearvars -except subjID biacID
close all;
end
