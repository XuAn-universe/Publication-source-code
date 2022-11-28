function result = analyze_performance(gttime, raise_interval, t)
%*---------------------------------------------------------------------*
% Huang Lab
% Duke University
% Author : Xu An, Jan 2022
% xu.an@duke.edu
% Version: 1.0
%*---------------------------------------------------------------------*
nmiss = 0;
nextra = 0;
nhit = 0;
terror = [];
for i = 1:size(raise_interval, 2)
    temp = find(gttime >= t(raise_interval(1, i)) & gttime < t(raise_interval(2, i)));
    if isempty(temp)
        nextra = nextra+1;
    else
        nhit = nhit+1;
        nmiss = nmiss+numel(temp)-1;
        terror_temp = t(raise_interval(3, i))-gttime(temp);
        [~, terrorID] = min(abs(terror_temp));
        terror = [terror; terror_temp(terrorID)];
        gttime(temp) = [];
    end
end
result.nmiss = nmiss+numel(gttime);
result.nextra = nextra;
result.nhit = nhit;
result.terror = terror;