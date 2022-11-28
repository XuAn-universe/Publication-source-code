%*---------------------------------------------------------------------*
% Huang Lab
% Duke University
% Author : Xu An, Feb 2022
% xu.an@duke.edu
% Version: 1.0
%*---------------------------------------------------------------------*
%%
% for HMM evaluation
curpwd = pwd;
try
   cd(pathname); 
end
[filename, pathname] = uigetfile('*.mat', 'Pick all of the data set', 'MultiSelect', 'on');
if isequal(filename, 0)
    cd(curpwd);
    return;
end

N = numel(filename);
hitrate_all = nan(1, N);
terror_each = cell(1, N);
terror_all = [];

for i = 1:N
    temp = load([pathname filename{i}]);
    result = temp.result;
    evalmatrix = result.evalmatrix;
    terror = result.terror;
    extrarate = evalmatrix(:, 1)./(evalmatrix(:, 1)+evalmatrix(:, 2));
    hitrate = evalmatrix(:, 2)./(evalmatrix(:, 3)+evalmatrix(:, 2));
    [hitrate_all(i), bestID] = max(hitrate);
    terror_each{i} = terror{bestID};
    terror_all = [terror_all; terror{bestID}];
end

figure;
bar(1, mean(hitrate_all), 1, 'FaceColor', [0.5 0.5 0.5]);
hold on;
plot(1, hitrate_all, 'o', 'MarkerFaceColor', [1 1 1], 'MarkerEdgeColor', [0 0 0]);
errorbar(1, mean(hitrate_all), std(hitrate_all)/sqrt(numel(hitrate_all)), 'k', 'LineStyle', 'none', 'CapSize', 15);
set(gca, 'XTick', 1, 'XTickLabel', '', 'XLim', [0 2], 'TickLength', [0 0], 'FontSize', 12);
ylabel('Hit rate');
box off;

[terror_min, terror_max] = bounds(terror_all);
tstep = 1/120*2;
edges = terror_min:tstep:terror_max;
count_all = nan(N, numel(edges)-1);
for i = 1:N
    count_all(i, :) = histcounts(terror_each{i}, edges, 'Normalization', 'probability');
end

figure;
bar(terror_min+tstep/2:tstep:terror_max-tstep/2, mean(count_all));
hold on;
errorbar(terror_min+tstep/2:tstep:terror_max-tstep/2, mean(count_all), [], std(count_all)/sqrt(size(count_all, 1)), 'k', 'LineStyle', 'none', 'CapSize', 4);
xlabel('Error (s)');
ylabel('Probability');
set(gca, 'TickDir', 'out', 'FontSize', 12);
box off;

figure;
histogram(terror_all, edges, 'Normalization', 'probability');
xlabel('Error (s)');
ylabel('Probability');
set(gca, 'TickDir', 'out', 'FontSize', 12);
box off;

prctile(abs(terror_all), 85)
terror_max = max(abs(terror_all));
edges = 0:tstep:terror_max;
count_all = nan(N, numel(edges)-1);
prctile_all = nan(1, N);
for i = 1:N
    count_all(i, :) = histcounts(abs(terror_each{i}), edges, 'Normalization', 'probability');
    prctile_all(i) = sum(abs(terror_each{i}) < 0.1)/numel(terror_each{i});
end
disp(['Error < 0.1 : ' num2str(mean(prctile_all)) '+-' num2str(std(prctile_all)/sqrt(N)) '; mean+-SEM']);
cumprob = cumsum(mean(count_all));
edges(find(cumprob >= 0.85, 1, 'first')+1)