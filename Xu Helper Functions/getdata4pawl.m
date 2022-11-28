function outputCell = getdata4pawl(inputCell)
%*---------------------------------------------------------------------*
% Huang Lab
% Duke University
% Author : Xu An, Jan 2022
% xu.an@duke.edu
% Version: 1.0
%*---------------------------------------------------------------------*
% function for DLtest
sig = inputCell{1};
sig(4, :) = [];
labels = inputCell{2};
labels(2, :) = [];
labels = removecats(labels);
outputCell = {sig, labels};
end