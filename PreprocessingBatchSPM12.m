clear;
spm fmri; % this is to open up a progress bar & print out graphics results

curDir = pwd;

load('dataAnalyzeBehavior','subjID','biacID'); % contains subject and scan IDs

for iS = 1 : length(subjID)
    scanInfo.subjID = subjID{iS}; 
    scanInfo.biacID = biacID{iS}; 
    scanInfo.dirData = fullfile(pwd,'Data'); % fMRI data directory
    scanInfo.dirSaveData = fullfile(scanInfo.dirData,scanInfo.subjID,scanInfo.biacID);
    scanInfo.dirAnat = fullfile(scanInfo.dirData,scanInfo.subjID,scanInfo.biacID,'Anat'); % location of T1 anatomical
    scanInfo.dirFunc = fullfile(scanInfo.dirData,scanInfo.subjID,scanInfo.biacID,'Func'); % location of functional scans
    scanInfo.expPrefix = 'SAT';
    scanInfo.TR = 2; % TR in sec.
    scanInfo.nTRDisdaq = 4; % number of TRs to discard
    scanInfo.nTRRun = [165 165 165 165 135 135 135 135 215 150]; % number of TRs for each session (not including disdaq)
    scanInfo.nRun = length(scanInfo.nTRRun); % total number of runs
    scanInfo.nSlice = 42; % number of slices (for slice timing correction)
    scanInfo.resolutionAnat = 1; % spatial resolution of T1 data
    scanInfo.resolutionFunc = 3; % spatial resolution of functional data
    scanInfo.smoothingKernelSize = 8; % size of Gaussian kernel for smoothing (e.g., 8)
    cd(scanInfo.dirSaveData);
    
    disp(['=========== Begin Preprocessing for ' scanInfo.subjID ' ==========']);
    
    % locate original functional files
    nameFunc = cell(1,scanInfo.nRun);
    prefixF = sprintf('^%s.*\\.img', scanInfo.expPrefix);
    for nR = 1 : scanInfo.nRun
        dirFuncRuns{nR} = fullfile(scanInfo.dirFunc,['Run' num2str(nR)]); % maps location of files for each run
        nameFiles{nR} = spm_select('List',dirFuncRuns{nR},prefixF); % detects .img file names
        for nT = 1 : scanInfo.nTRRun(nR)
            tempName = fullfile(dirFuncRuns{nR},nameFiles{nR}(nT+scanInfo.nTRDisdaq,:));
            nameFunc{1,nR}(nT,:) = {[tempName ',1']};
        end
    end
    
    % loacate T1 anantomical file 
    prefixA = sprintf('^%s.*\\.img', scanInfo.expPrefix);
    nameFileAnat = spm_select('List',scanInfo.dirAnat,prefixA);
    nameAnat = {[fullfile(scanInfo.dirAnat,nameFileAnat) ',1']}; % maps full directory of T1 image
    
    %% REALIGNMENT (ESTIMATE & RESLICE)
    % Motion correction
    % Do realignment and then slice timing
    % see http://imaging.mrc-cbu.cam.ac.uk/imaging/SliceTimingRikEmail for a discussion on this issue.
    
    matlabbatch{1}.spm.spatial.realign.estwrite.data = nameFunc;
    matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.quality = 0.9;
    matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.sep = 4;
    matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.fwhm = 5;
    matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.rtm = 1; % register to mean (don't change it!)
    matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.interp = 2;
    matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.wrap = [0 0 0];
    matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.weight = '';
    
    matlabbatch{1}.spm.spatial.realign.estwrite.roptions.which = [0 1]; % mean images only
    matlabbatch{1}.spm.spatial.realign.estwrite.roptions.interp = 4;
    matlabbatch{1}.spm.spatial.realign.estwrite.roptions.wrap = [0 0 0];
    matlabbatch{1}.spm.spatial.realign.estwrite.roptions.mask = 0; % no masking
    matlabbatch{1}.spm.spatial.realign.estwrite.roptions.prefix = 'r';
    
    %% SLICE TIMING CORRECTION
    
    matlabbatch{2}.spm.temporal.st.scans = nameFunc;
    matlabbatch{2}.spm.temporal.st.nslices = scanInfo.nSlice;
    matlabbatch{2}.spm.temporal.st.tr = scanInfo.TR;
    matlabbatch{2}.spm.temporal.st.ta = scanInfo.TR - (scanInfo.TR/scanInfo.nSlice);
    matlabbatch{2}.spm.temporal.st.so = [1:2:scanInfo.nSlice, 2:2:scanInfo.nSlice]; % slice order for Slice Timing
    matlabbatch{2}.spm.temporal.st.refslice = 1;
    matlabbatch{2}.spm.temporal.st.prefix = 'a';
    
    %% COREGISTRATION (ESTIMATE)
    % Overlay structural and functional images: Link functional scans to anatomical scan
    % REFERENCE IMAGE: read in mean functional image (remains stationary)
    % SOURCE IMAGE: read in anatomical T1 image (image to be jiggled)
    
    % locate mean functional file (reference image)
    nameFuncMean = {[fullfile(scanInfo.dirFunc,'Run1',['mean' nameFiles{1}(1+scanInfo.nTRDisdaq,:)]) ',1']};
    
    matlabbatch{3}.spm.spatial.coreg.estimate.ref = nameFuncMean;
    matlabbatch{3}.spm.spatial.coreg.estimate.source = nameAnat;
    matlabbatch{3}.spm.spatial.coreg.estimate.other = {''};
    matlabbatch{3}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi'; % Normalised Mutual Information
    matlabbatch{3}.spm.spatial.coreg.estimate.eoptions.sep = [4 2];
    matlabbatch{3}.spm.spatial.coreg.estimate.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
    matlabbatch{3}.spm.spatial.coreg.estimate.eoptions.fwhm = [7 7];
    
    %% SEGMENTATION
    % SPM12 uses deformation field produced in this step for normalization
    
    matlabbatch{4}.spm.spatial.preproc.channel.vols = nameAnat;
    matlabbatch{4}.spm.spatial.preproc.channel.biasreg = 0.001;
    matlabbatch{4}.spm.spatial.preproc.channel.biasfwhm = 60;
    matlabbatch{4}.spm.spatial.preproc.channel.write = [0 1]; % save bias corrected
    matlabbatch{4}.spm.spatial.preproc.tissue(1).tpm = {'C:\Users\ho20\Documents\MATLAB\SPM12\tpm\TPM.nii,1'};
    matlabbatch{4}.spm.spatial.preproc.tissue(1).ngaus = 1;
    matlabbatch{4}.spm.spatial.preproc.tissue(1).native = [1 0];
    matlabbatch{4}.spm.spatial.preproc.tissue(1).warped = [0 0];
    matlabbatch{4}.spm.spatial.preproc.tissue(2).tpm = {'C:\Users\ho20\Documents\MATLAB\SPM12\tpm\TPM.nii,2'};
    matlabbatch{4}.spm.spatial.preproc.tissue(2).ngaus = 1;
    matlabbatch{4}.spm.spatial.preproc.tissue(2).native = [1 0];
    matlabbatch{4}.spm.spatial.preproc.tissue(2).warped = [0 0];
    matlabbatch{4}.spm.spatial.preproc.tissue(3).tpm = {'C:\Users\ho20\Documents\MATLAB\SPM12\tpm\TPM.nii,3'};
    matlabbatch{4}.spm.spatial.preproc.tissue(3).ngaus = 2;
    matlabbatch{4}.spm.spatial.preproc.tissue(3).native = [1 0];
    matlabbatch{4}.spm.spatial.preproc.tissue(3).warped = [0 0];
    matlabbatch{4}.spm.spatial.preproc.tissue(4).tpm = {'C:\Users\ho20\Documents\MATLAB\SPM12\tpm\TPM.nii,4'};
    matlabbatch{4}.spm.spatial.preproc.tissue(4).ngaus = 3;
    matlabbatch{4}.spm.spatial.preproc.tissue(4).native = [1 0];
    matlabbatch{4}.spm.spatial.preproc.tissue(4).warped = [0 0];
    matlabbatch{4}.spm.spatial.preproc.tissue(5).tpm = {'C:\Users\ho20\Documents\MATLAB\SPM12\tpm\TPM.nii,5'};
    matlabbatch{4}.spm.spatial.preproc.tissue(5).ngaus = 4;
    matlabbatch{4}.spm.spatial.preproc.tissue(5).native = [1 0];
    matlabbatch{4}.spm.spatial.preproc.tissue(5).warped = [0 0];
    matlabbatch{4}.spm.spatial.preproc.tissue(6).tpm = {'C:\Users\ho20\Documents\MATLAB\SPM12\tpm\TPM.nii,6'};
    matlabbatch{4}.spm.spatial.preproc.tissue(6).ngaus = 2;
    matlabbatch{4}.spm.spatial.preproc.tissue(6).native = [0 0];
    matlabbatch{4}.spm.spatial.preproc.tissue(6).warped = [0 0];
    matlabbatch{4}.spm.spatial.preproc.warp.mrf = 1;
    matlabbatch{4}.spm.spatial.preproc.warp.cleanup = 1;
    matlabbatch{4}.spm.spatial.preproc.warp.reg = [0 0.001 0.5 0.05 0.2];
    matlabbatch{4}.spm.spatial.preproc.warp.affreg = 'mni';
    matlabbatch{4}.spm.spatial.preproc.warp.fwhm = 0;
    matlabbatch{4}.spm.spatial.preproc.warp.samp = 3;
    matlabbatch{4}.spm.spatial.preproc.warp.write = [0 1]; % FORWARD: for spatially normalizing images to mni space, need to use the forward deformation
    
    %% NORMALIZE (WRITE)
    % Warp images to fit to a standard template brain (MNI)
    
    % locate deformation field (created as y*.nii)
    nameFileDeformField = ['y_' nameFileAnat(1:end-4) '.nii'];
    nameDeformField = {fullfile(scanInfo.dirAnat,nameFileDeformField)}; 
    
    % locate slice time corrected functional images
    nImage = 1;
    for nR = 1 : scanInfo.nRun
        for nT = 1 : scanInfo.nTRRun(nR)
            nameFuncA(nImage,:) = {[fullfile(dirFuncRuns{nR},['a' nameFiles{nR}(nT+scanInfo.nTRDisdaq,:)]) ',1']};
            nImage = nImage + 1;
        end
    end
    
    matlabbatch{5}.spm.spatial.normalise.write.subj.def = nameDeformField;
    matlabbatch{5}.spm.spatial.normalise.write.subj.resample = nameFuncA;
    matlabbatch{5}.spm.spatial.normalise.write.woptions.bb = [-78 -112 -70; 78 76 85];
    matlabbatch{5}.spm.spatial.normalise.write.woptions.vox = repmat(scanInfo.resolutionFunc,1,3); % resample as the original resolution
    matlabbatch{5}.spm.spatial.normalise.write.woptions.interp = 4;
    matlabbatch{5}.spm.spatial.normalise.write.woptions.prefix = 'w';
    
    % T1 Anatomical Rewrite
    % This is to view the quality of normalization using T1 (better resolution than functional scans)    
    % locate bias corrected T1 structural image
    nameFileAnatCorrected = ['m' nameFileAnat(1:end-4) '.nii'];
    nameAnatCorrected = {[fullfile(scanInfo.dirAnat,nameFileAnatCorrected) ',1']}; % maps full directory of T1 image
   
    matlabbatch{6}.spm.spatial.normalise.write.subj.def = nameDeformField;
    matlabbatch{6}.spm.spatial.normalise.write.subj.resample = nameAnatCorrected;
    matlabbatch{6}.spm.spatial.normalise.write.woptions.bb = [-78 -112 -70; 78 76 85];
    matlabbatch{6}.spm.spatial.normalise.write.woptions.vox = repmat(scanInfo.resolutionAnat,1,3); % resample as the original resolution
    matlabbatch{6}.spm.spatial.normalise.write.woptions.interp = 4;
    matlabbatch{6}.spm.spatial.normalise.write.woptions.prefix = 'w';

    %% SMOOTHING
    
    % locate normalized functional images
    nImage = 1;
    for nR = 1 : scanInfo.nRun
        for nT = 1 : scanInfo.nTRRun(nR)
            nameFuncW(nImage,:) = {[fullfile(dirFuncRuns{nR},['wa' nameFiles{nR}(nT+scanInfo.nTRDisdaq,:)]) ',1']};
            nImage = nImage + 1;
        end
    end
    
    matlabbatch{7}.spm.spatial.smooth.data = nameFuncW;
    matlabbatch{7}.spm.spatial.smooth.fwhm = repmat(scanInfo.smoothingKernelSize,1,3);
    matlabbatch{7}.spm.spatial.smooth.dtype = 0;
    matlabbatch{7}.spm.spatial.smooth.im = 0;
    matlabbatch{7}.spm.spatial.smooth.prefix = 's';
    
    %% RUN SPM
    
    spm('defaults', 'FMRI');
    spm_jobman('serial',matlabbatch)
    
    %% DONE
    
    save([scanInfo.subjID '_paramPreprocess.mat'],'scanInfo','nameFiles','nameFunc','nameFileAnat','nameAnat');
    
    % copy and rear
    
    disp(['=========== Done Preprocessing: ' scanInfo.subjID ' ==========']);
    
    clearvars -EXCEPT subjID biacID curDir    
end

cd(curDir);
close all;
