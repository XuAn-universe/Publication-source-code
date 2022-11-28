function result = time2category(t, trange, categorynames)
%*---------------------------------------------------------------------*
% Huang Lab
% Duke University
% Author : Xu An, Jan 2022
% xu.an@duke.edu
% Version: 1.0
%*---------------------------------------------------------------------*
result = zeros(1, numel(t));
for i = 1:numel(trange)
    tpair = trange{i};
    for j = 1:size(tpair, 1)
        result(t >= tpair(j, 1) & t <= tpair(j, 2)) = i;
    end
end
result = categorical(result, 0:numel(categorynames), [{'Unlabelled'} categorynames]);