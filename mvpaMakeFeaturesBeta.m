function features = mvpaMakeFeaturesBeta(nameIdx,dirData,betaID,saveName)
% this function extracts SPM beta maps (parameter estimates)
% nameIdx: full path to the neighbor file (searchlight sphere file)
% dirData: GLM output file directory 
% betaID: T map number, e.g., SPMT_0001.nii 
% saveName: output file name 
% flagMultipleGLM: input number of GLM analysis 

% load neighbor data (e.g., searchIdx2.dat - searchlight with radius = 2 voxels)
fid = fopen(nameIdx, 'r');
len = fread(fid, 1, 'int32');
sL2Voxel = fread(fid, len, 'int32'); %grey matter voxel (because search*.dat files are specified for 53x63x46 voxels, need to process each subject's data with this specification; use rwa*.img files)

% load SPM beta maps (parameter estimates)
features = [];
for iR = 1 : size(betaID,1)
    for iC = betaID(iR,:)
        fileName = fullfile(dirData,sprintf('beta_%04d.nii', iC)); % load parameter estimates, beta
        disp(fileName);
        fidB = fopen(fileName, 'r');
        data = fread(fidB, inf, 'float32');
        fclose(fidB);
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

if exist('saveName','var')
    save(saveName, 'features'); % save preprocessed images for searchlight analysis
end

fclose(fid);