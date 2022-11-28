%*---------------------------------------------------------------------*
% Huang Lab
% Duke University
% Author : Xu An, Oct 2021
% xu.an@duke.edu
% Version: 1.0
%*---------------------------------------------------------------------*
%% Fisher's exact test, chisquare test is crosstab
x = table([208; 105], [0; 42], 'VariableNames', {'success', 'failure'}, 'RowNames', {'saline', 'muscimol'})
[h, p, stats] = fishertest(x, 'Tail', 'both', 'Alpha', 0.05)