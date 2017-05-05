clear;
addpath('C:\Users\ho20\Desktop\fMRI_Data\SAT.01\dataBehavioral');
load('dataExtractBehavior3','scanP');
dirData = 'C:\Users\ho20\Desktop\fMRI_Data\SAT.01\dataFMRI';
curDir = pwd;

nameAnalysis = 'MVPA_CueValue_pmodFeedbackRT_perRun_NativeSpace';
nameHdrPrefix = {'cvNP_CueRegression_','cvTP_CueRegression_','mseNP_CueRegression_','mseTP_CueRegression_'};
foldername = 'C:\Users\ho20\Desktop\fMRI_Data\SAT.01\MVPA\AccuracyMap\MVPA_CueValueRegression_WithinSubject_NativeSpace\Tmap_searchIdx4';

for iS = 1 : scanP.nSubj
    dirAnat = fullfile(dirData,scanP.subjID{iS},scanP.biacID{iS},'Anat'); % T1 anatomical
    dirFunc = fullfile(dirData,scanP.subjID{iS},scanP.biacID{iS},'Func'); % functional scans
    
    % resample the grey matter mask to match functional space
%     nameGM = fullfile(dirAnat,['c1' scanP.subjID{iS} '_anat0001.nii']);
%     nameTargetSpace = fullfile(dirFunc,'Run1',['ar' scanP.subjID{iS} '_Run1_0005.img']);    
%     nameImages = {nameTargetSpace, nameGM};
%     optReslice.which = [1 0]; % images 2...n (don't reslice the sample image)
%     optReslice.interp = 0; % nearest neighbor
%     optReslice.wrap = [0 0 0];
%     optReslice.mask = 1;
%     optReslice.prefix = 'r';    
%     spm_reslice(nameImages,optReslice);
    
    % create searchlight index file
    nameGMResampled = fullfile(dirAnat,['frc1' scanP.subjID{iS} '_anat0001.img']);
    dim =[64 64 42]; 
    threshold = 1;
    r = 5; 
    saveName = fullfile(curDir,'NativeSpace',[scanP.subjID{iS} '_searchIdx' num2str(r) '.dat']); 
    GenerateIdxFile(nameGMResampled, dim, saveName, threshold, r);

    % copy spmT files to create header files
%     dirDataT = fullfile(dirData,scanP.subjID{iS},'Analysis',nameAnalysis,'spmT_0001.nii');
%     copyfile(dirDataT,fullfile(curDir,'NativeSpace',[scanP.subjID{iS} '_spmT_0001.nii']));   

%     nameHdr = fullfile(curDir,'NativeSpace',['f' scanP.subjID{iS} '_spmT_0001.hdr']);
%     for iH = 1 : 4
%         copyfile(nameHdr,fullfile(foldername,[nameHdrPrefix{iH} scanP.subjID{iS} '.hdr']));  
%     end
end
