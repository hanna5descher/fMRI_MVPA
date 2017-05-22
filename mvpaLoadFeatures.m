function [f1, f2] = mvpaLoadFeatures(dirFeatures,subjID,groupLabel1,groupLabel2)
% dirFeatures: directory of feature files
% subjID: list of subject IDs to be included
% groupLabel1 & 2: row index of features to be loaded

f1 = []; f2 = [];
for i = 1 : length(subjID)
    tempF = load(fullfile(dirFeatures,subjID{i}),'features');
    for nF = 1 : length(groupLabel1)
        f1 = [f1; tempF.features(groupLabel1(nF),:)];
        f2 = [f2; tempF.features(groupLabel2(nF),:)];
    end    
    clear tempF
end