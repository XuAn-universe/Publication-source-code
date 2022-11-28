function matchedIDs = matchlabels(labels, sparsename)
%*---------------------------------------------------------------------*
% Huang Lab
% Duke University
% Author : Xu An, Jan 2022
% xu.an@duke.edu
% Version: 1.0
%*---------------------------------------------------------------------*
% function to match the number of samples for imbalanced classes
folds = 4; % unlabelled/labelled
temp = labels == sparsename;
temp = [0 temp 0];
startID = find(diff(temp) == 1);
endID = find(diff(temp) == -1)-1;
endID = endID+(endID-startID+1)*folds;
endID(endID > numel(labels)) = numel(labels);
matchedIDs = false(1, numel(labels));
for i = 1:numel(startID)
    matchedIDs(startID(i):endID(i)) = true;
end