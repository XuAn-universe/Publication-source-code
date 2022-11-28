%*---------------------------------------------------------------------*
% Huang Lab
% Duke University
% Author : Xu An, Feb 2022
% xu.an@duke.edu
% Version: 1.0
%*---------------------------------------------------------------------*
%%
fezf1 = fezf1/(sum(fezf1(:)));
fezf2 = fezf2/(sum(fezf2(:)));
fezf3 = fezf3/(sum(fezf3(:)));
fezf(:, 1) = mean([fezf1(:, 1) fezf2(:, 1) fezf3(:, 1)], 2);
fezf(:, 2) = mean([fezf1(:, 2) fezf2(:, 2) fezf3(:, 2)], 2);

%%
plexin1 = plexin1/(sum(plexin1(:)));
plexin2 = plexin2/(sum(plexin2(:)));
plexin(:, 1) = mean([plexin1(:, 1) plexin2(:, 1)], 2);
plexin(:, 2) = mean([plexin1(:, 2) plexin2(:, 2)], 2);