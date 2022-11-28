%*---------------------------------------------------------------------*
% Huang Lab
% Duke University
% Author : Xu An, May 2022
% xu.an@duke.edu
% Version: 1.0
%*---------------------------------------------------------------------*
%%
% combine result from different sessions
result.adjustment2firstbite = [result1.adjustment2firstbite; result2.adjustment2firstbite];
result.retrieval2sit = [result1.retrieval2sit; result2.retrieval2sit];
result.retrieval2lefthandreach = [result1.retrieval2lefthandreach; result2.retrieval2lefthandreach];
result.retrieval2righthandreach = [result1.retrieval2righthandreach; result2.retrieval2righthandreach];

%%
% combine result from different sessions
result.adjustment2firstbite = [result1.adjustment2firstbite; result2.adjustment2firstbite];
result.retrieval2sit = [result1.retrieval2sit; result2.retrieval2sit];
result.retrieval2lefthandreach = [result1.retrieval2lefthandreach; result2.retrieval2lefthandreach];
result.retrieval2righthandreach = [result1.retrieval2righthandreach; result2.retrieval2righthandreach];
result.firstbite2withdraw = [result1.firstbite2withdraw; result2.firstbite2withdraw];
result.firstadjustment2withdraw = [result1.firstadjustment2withdraw; result2.firstadjustment2withdraw];

%%
% combine result from different sessions
names = fieldnames(result1);
for i = 1:numel(names)
    eval(['result.' names{i} ' = [result1.' names{i} '; result2.' names{i} '];']);
end

%% feeding duration
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

feedingduration_all = [];
feedingduration_each = nan(1, N);

for i = 1:N
    temp = load([pathname filename{i}]);
    result = temp.result;
    feedingduration_all = [feedingduration_all; result.feedingduration];
    feedingduration_each(i) = mean(result.feedingduration);
end

cd(curpwd);

%%
figure;
paired_plot(feedingduration_each_saline, feedingduration_each_muscimol, 'Feeding duration (s)', filename);

%% adjustment2firstbite, all2retrievalend
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

tstep = 0.1;
tstep_s = 1/120;

adjustment2firstbite_all = [];
adjustment2firstbite_each = cell(1, N);
propnegtive = nan(1, N);
sit2retrievalend_all = [];
sit2retrievalend_each = cell(1, N);
lefthandreach2retrievalend_all = [];
lefthandreach2retrievalend_each = cell(1, N);
righthandreach2retrievalend_all = [];
righthandreach2retrievalend_each = cell(1, N);

for i = 1:N
    temp = load([pathname filename{i}]);
    result = temp.result;
    adjustment2firstbite_all = [adjustment2firstbite_all; result.adjustment2firstbite];
    adjustment2firstbite_each{i} = result.adjustment2firstbite;
    propnegtive(i) = sum(result.adjustment2firstbite < 0)/numel(result.adjustment2firstbite);
    sit2retrievalend_all = [sit2retrievalend_all; result.sit2retrievalend];
    sit2retrievalend_each{i} = result.sit2retrievalend;
    lefthandreach2retrievalend_all = [lefthandreach2retrievalend_all; result.lefthandreach2retrievalend];
    lefthandreach2retrievalend_each{i} = result.lefthandreach2retrievalend;
    righthandreach2retrievalend_all = [righthandreach2retrievalend_all; result.righthandreach2retrievalend];
    righthandreach2retrievalend_each{i} = result.righthandreach2retrievalend;
end
disp(['adjustment2firstbite < 0 : ' num2str(mean(propnegtive)) '+-' num2str(std(propnegtive)/sqrt(N)) '; mean+-SEM']);

make_histograms(adjustment2firstbite_all, tstep, adjustment2firstbite_each, 'Time (s)', 'adjustment2firstbite');
make_histograms(sit2retrievalend_all, 2*tstep_s, sit2retrievalend_each, 'Time (s)', 'sit2retrievalend');
make_histograms(lefthandreach2retrievalend_all, 2*tstep_s, lefthandreach2retrievalend_each, 'Time (s)', 'lefthandreach2retrievalend');
make_histograms(righthandreach2retrievalend_all, 2*tstep_s, righthandreach2retrievalend_each, 'Time (s)', 'righthandreach2retrievalend');

cd(curpwd);

%% all2retrievalstart
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

tstep = 0.1;
tstep_s = 1/120;

sit2retrievalstart_all = [];
sit2retrievalstart_each = cell(1, N);
lefthandreach2retrievalstart_all = [];
lefthandreach2retrievalstart_each = cell(1, N);
righthandreach2retrievalstart_all = [];
righthandreach2retrievalstart_each = cell(1, N);

for i = 1:N
    temp = load([pathname filename{i}]);
    result = temp.result;
    sit2retrievalstart_all = [sit2retrievalstart_all; result.sit2retrievalstart];
    sit2retrievalstart_each{i} = result.sit2retrievalstart;
    lefthandreach2retrievalstart_all = [lefthandreach2retrievalstart_all; result.lefthandreach2retrievalstart];
    lefthandreach2retrievalstart_each{i} = result.lefthandreach2retrievalstart;
    righthandreach2retrievalstart_all = [righthandreach2retrievalstart_all; result.righthandreach2retrievalstart];
    righthandreach2retrievalstart_each{i} = result.righthandreach2retrievalstart;
end

make_histograms(sit2retrievalstart_all, 2*tstep_s, sit2retrievalstart_each, 'Time (s)', 'sit2retrievalstart', [0 0.35]);
make_histograms(lefthandreach2retrievalstart_all, 2*tstep_s, lefthandreach2retrievalstart_each, 'Time (s)', 'lefthandreach2retrievalstart', [0 0.35]);
make_histograms(righthandreach2retrievalstart_all, 2*tstep_s, righthandreach2retrievalstart_each, 'Time (s)', 'righthandreach2retrievalstart', [0 0.35]);

cd(curpwd);

%% all2withdraw
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

tstep = 0.05;
tstep_s = 1/120;

firstadjustment2withdraw_all = [];
firstadjustment2withdraw_each = cell(1, N);
firstbite2withdraw_all = [];
firstbite2withdraw_each = cell(1, N);

for i = 1:N
    temp = load([pathname filename{i}]);
    result = temp.result;
    firstadjustment2withdraw_all = [firstadjustment2withdraw_all; result.firstadjustment2withdraw];
    firstadjustment2withdraw_each{i} = result.firstadjustment2withdraw;
    firstbite2withdraw_all = [firstbite2withdraw_all; result.firstbite2withdraw];
    firstbite2withdraw_each{i} = result.firstbite2withdraw;
end

firstadjustment2withdraw_count_all = make_histograms(firstadjustment2withdraw_all, tstep, firstadjustment2withdraw_each, 'Time (s)', 'firstadjustment2withdraw', [0 1]);
firstbite2withdraw_count_all = make_histograms(firstbite2withdraw_all, tstep, firstbite2withdraw_each, 'Time (s)', 'firstbite2withdraw', [0 1]);
figure;
hp(1) = plot_tj_MeanSEM(0+tstep/2:tstep:1-tstep/2, firstadjustment2withdraw_count_all', [1 0 0], [1 0 0], 'Time (s)', 'Probability', '');
hp(2) = plot_tj_MeanSEM(0+tstep/2:tstep:1-tstep/2, firstbite2withdraw_count_all', [0 1 0], [0 1 0], 'Time (s)', 'Probability', '');
xlim([0 1]);
legend(hp, {'firstadjustment2withdraw', 'firstbite2withdraw'}, 'Location', 'NorthEast');
legend('boxoff');

cd(curpwd);

%% inter bout interval
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

tstep = 1/120;

interbout_interval_all = [];
interbout_interval_each = cell(1, N);

for i = 1:N
    temp = load([pathname filename{i}]);
    result = temp.result;
    interbout_interval_all = [interbout_interval_all; result.interbout_interval];
    interbout_interval_each{i} = result.interbout_interval;
end

make_histograms(interbout_interval_all, tstep, interbout_interval_each, 'Time (s)', 'Inter bout interval');

cd(curpwd);

function count_all = make_histograms(data_all, tstep, data_each, xlabel_text, title_text, tbounds)
if nargin < 6
    tbounds = [NaN NaN];
end
N = numel(data_each);
[t_min, t_max] = bounds(data_all);
if ~isnan(tbounds(1))
    t_min = tbounds(1);
end
if ~isnan(tbounds(2))
    t_max = tbounds(2);
end
edges = t_min:tstep:t_max;
count_all = nan(N, numel(edges)-1);
for i = 1:N
    count_all(i, :) = histcounts(data_each{i}, edges, 'Normalization', 'probability');
end

figure;
bar(t_min+tstep/2:tstep:t_max-tstep/2, mean(count_all));
hold on;
errorbar(t_min+tstep/2:tstep:t_max-tstep/2, mean(count_all), [], std(count_all)/sqrt(size(count_all, 1)), 'k', 'LineStyle', 'none', 'CapSize', 0);
xlabel(xlabel_text);
ylabel('Probability');
title(title_text);
set(gca, 'TickDir', 'out', 'FontSize', 12);
box off;

figure;
histogram(data_all, edges, 'Normalization', 'probability');
xlabel(xlabel_text);
ylabel('Probability');
title(title_text);
set(gca, 'TickDir', 'out', 'FontSize', 12);
box off;
end