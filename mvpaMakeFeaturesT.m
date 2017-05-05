function features = mvpaMakeFeaturesT(subjID,nameIdx,dirData,contrastID,saveName,flagMultipleGLM)
% subjID: subject ID, e.g., SAT01
% nameIdx: full path to the neighbor file (searchlight sphere file)
% dirData: GLM output file directory
% contrastID: T map number, e.g., SPMT_0001.nii 
% saveName: output file name 
% flagMultipleGLM: input number of GLM analysis 

% this function preprocesses SPM T maps for MVPA
if nargin < 6
    flagMultipleGLM = 1;
    dirCont{1} = dirData;
else
    dirCont = dirData;
end

% load searchlight data 
fid = fopen(nameIdx, 'r');
len = fread(fid, 1, 'int32');
sL2Voxel = fread(fid, len, 'int32'); 

% load SPM T maps (contrasts)
features = [];
for iC = contrastID
    for iG = 1 : flagMultipleGLM
        fileName = fullfile(dirCont{iG},sprintf('spmT_%04d.nii', iC)); % load T map (image dimension must match the dimension of the searchlight file)
        disp(fileName);
        fidT = fopen(fileName, 'r');
        data = fread(fidT, inf, 'float32');
        fclose(fidT);
        data = data(89:end); % extract data that are equivalent to .img format
        
        features = [features; data(sL2Voxel)']; % read in T maps and reorganize them to selected grey matter voxels
    end
end

disp(['Feature matrix size: ' num2str(size(features,1)) ' by ' num2str(size(features,2))]);

%center/normalize data for each voxel
for iF = 1 : size(features, 2)
    x = features(:,iF);
    if sum(abs(x)) > 0
        x = (x - mean(x)) / std(x);
        %x = (x - min(x)) / (max(x) - min(x));
        %x = (x - mean(x));
        features(:,iF) = x(:);
    end
end

%save(saveName, 'features'); % save preprocessed images for searchlight analysis
fclose(fid); 
