%*---------------------------------------------------------------------*
% Huang Lab
% Duke University
% Author : Xu An, Feb 2022
% xu.an@duke.edu
% Version: 1.0
%*---------------------------------------------------------------------*
%%
% combine result from different sessions
result = [];
names = fieldnames(result1);
for i = 1:numel(names)
    if length(names{i}) > 6 && strcmp(names{i}(1:6), 'zscore')
        eval(['result.' names{i} ' = [result1.' names{i} ' result2.' names{i} '];']);
    elseif ~strcmp(names{i}, 't')
        eval(['result.' names{i} ' = [result1.' names{i} '; result2.' names{i} '];']);
    end
end
result.t = result1.t;

%%
% combine result from different sessions
result.zscore_retrievalstart = [result1.zscore_retrievalstart result2.zscore_retrievalstart];
result.zscore_biteboutstart = [result1.zscore_biteboutstart result2.zscore_biteboutstart];
result.zscore_biteboutstartAdj = [result1.zscore_biteboutstartAdj result2.zscore_biteboutstartAdj];
result.zscore_biteboutstartNAdj = [result1.zscore_biteboutstartNAdj result2.zscore_biteboutstartNAdj];
result.zscore_barriercross = [result1.zscore_barriercross result2.zscore_barriercross];
result.zscore_lastbite = [result1.zscore_lastbite result2.zscore_lastbite];
result.t = result1.t;

%%
% combine result from different sessions for lick
result.zscore_lickstart = [result1.zscore_lickstart result2.zscore_lickstart];
result.tlick = result1.t;

%%
% combine result from different sessions for bite device
result.zscore_bite = [result1.zscore_bite result2.zscore_bite];
result.zscore_firstbite = [result1.zscore_firstbite result2.zscore_firstbite];
result.tbite = result1.t;

%%
% compute peak to trough difference directly from plots
p.Position(2)-t.Position(2)

%%
% for 15 mm pasta feeding analysis (old)
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

nchannel = 2;
trange_retrieval = [-1 1];
colors = [0.3010 0.7450 0.9330; 0.8500 0.3250 0.0980];
hp = zeros(1, nchannel);
for i = 1:nchannel
    lgdtext{i} = ['Channel ' num2str(i)];
end

dzscore_retrieval = nan(N, nchannel);
dzscore_biteboutstart = nan(N, nchannel);
peakzscore_biteboutstart = nan(N, nchannel);

for i = 1:N
    temp = load([pathname filename{i}]);
    result = temp.result;
    t = result.t;
    zscore_retrievalstart = result.zscore_retrievalstart;
    figure('Name', filename{i});
    subplot(1, 3, 1);
    for j = 1:nchannel
        hp(j) = plot_tj_MeanSEM(t, zscore_retrievalstart(:, :, j), colors(j, :), colors(j, :), 'Time (s)', 'Z-score', 'aligned to Retrieval Start with Mouth');
    end
    yl = ylim;
    plot([0 0], yl, '--k', 'LineWidth', 1);
    legend(hp, lgdtext, 'Location', 'NorthEast');
    legend('boxoff');
    for j = 1:nchannel
        zscore_avg = mean(zscore_retrievalstart(:, :, j), 2, 'omitnan');
        [zscore_max, maxID] = max(zscore_avg(t > 0 & t < trange_retrieval(2)));
        temp = t(t > 0 & t < trange_retrieval(2));
        dzscore_retrieval(i, j) = zscore_max-min(zscore_avg(t < temp(maxID) & t > trange_retrieval(1)));
    end
    
    zscore_biteboutstart = result.zscore_biteboutstart;
    subplot(1, 3, 2);
    for j = 1:nchannel
        hp(j) = plot_tj_MeanSEM(t, zscore_biteboutstart(:, :, j), colors(j, :), colors(j, :), 'Time (s)', 'Z-score', 'aligned to Bite Bout Start');
    end
    yl = ylim;
    plot([0 0], yl, '--k', 'LineWidth', 1);
    legend(hp, lgdtext, 'Location', 'NorthEast');
    legend('boxoff');
    for j = 1:nchannel
        zscore_avg = mean(zscore_biteboutstart(:, :, j), 2, 'omitnan');
        [zscore_min, minID] = min(zscore_avg(t >= 0 & t < trange_retrieval(2)));
        temp = t(t >= 0 & t < trange_retrieval(2));
        dzscore_biteboutstart(i, j) = max(zscore_avg(t > temp(minID)))-zscore_min;
        
        peakzscore_biteboutstart(i, j) = max(zscore_avg(t > 0));
    end
    
    zscore_lastbite = result.zscore_lastbite;
    subplot(1, 3, 3);
    for j = 1:nchannel
        hp(j) = plot_tj_MeanSEM(t, zscore_lastbite(:, :, j), colors(j, :), colors(j, :), 'Time (s)', 'Z-score', 'aligned to Chew Start');
    end
    yl = ylim;
    plot([0 0], yl, '--k', 'LineWidth', 1);
    legend(hp, lgdtext, 'Location', 'NorthEast');
    legend('boxoff');
end
figure('Name', 'Retrieval Start');
disp('\Delta Z-score retrieval');
paired_plot(dzscore_retrieval(:, 1), dzscore_retrieval(:, 2), '\Delta Z-score', filename, 'bar');
set(gca, 'XTick', [1 2], 'XTickLabel', {'CFA', 'RFO'});

figure('Name', 'Bite Bout Start');
disp('\Delta Z-score biteboutstart');
paired_plot(dzscore_biteboutstart(:, 1), dzscore_biteboutstart(:, 2), '\Delta Z-score', filename, 'dot');
set(gca, 'XTick', [1 2], 'XTickLabel', {'CFA', 'RFO'});

figure('Name', 'Bite Bout Start');
disp('Peak Z-score biteboutstart');
paired_plot(peakzscore_biteboutstart(:, 1), peakzscore_biteboutstart(:, 2), 'Peak Z-score', filename, 'dot');
set(gca, 'XTick', [1 2], 'XTickLabel', {'CFA', 'RFO'});

cd(curpwd);

%%
% for 15 mm pasta feeding analysis (new)
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

nchannel = 2;
trange_retrieval = [-1 1];
colors = [0.3010 0.7450 0.9330; 0.8500 0.3250 0.0980];
hp = zeros(1, nchannel);
for i = 1:nchannel
    lgdtext{i} = ['Channel ' num2str(i)];
end

zscore_avg_handadj = [];
zscore_avg_bite = [];
zscore_avg_chew = [];

dzscore_handadj = nan(N, nchannel);
dzscore_bite = nan(N, nchannel);
dzscore_chew = nan(N, nchannel);
peakzscore_handadj = nan(N, nchannel);
peakzscore_bite = nan(N, nchannel);
peakzscore_chew = nan(N, nchannel);

for i = 1:N
    temp = load([pathname filename{i}]);
    result = temp.result;
    t = result.t;
    zscore_handadj = result.zscore_handadj;
    zscore_bite = result.zscore_bite;
    zscore_chew = result.zscore_chew;
    figure('Name', filename{i});
    subplot(1, 3, 1);
    for j = 1:nchannel
        hp(j) = plot_tj_MeanSEM(t, zscore_handadj(:, :, j), colors(j, :), colors(j, :), 'Time (s)', 'Z-score', 'aligned to Hand Adjustment');
    end
    yl = ylim;
    plot([0 0], yl, '--k', 'LineWidth', 1);
    legend(hp, lgdtext, 'Location', 'NorthEast');
    legend('boxoff');
    for j = 1:nchannel
        zscore_avg = mean(zscore_handadj(:, :, j), 2, 'omitnan');
        zscore_avg_handadj(:, i, j) = zscore_avg;
%         [zscore_max, maxID] = max(zscore_avg(t > 0 & t < trange_retrieval(2)));
%         temp = t(t > 0 & t < trange_retrieval(2));
%         dzscore_handadj(i, j) = zscore_max-min(zscore_avg(t < temp(maxID) & t > trange_retrieval(1)));
    end
    
    zscore_bite = result.zscore_bite;
    subplot(1, 3, 2);
    for j = 1:nchannel
        hp(j) = plot_tj_MeanSEM(t, zscore_bite(:, :, j), colors(j, :), colors(j, :), 'Time (s)', 'Z-score', 'aligned to Bite');
    end
    yl = ylim;
    plot([0 0], yl, '--k', 'LineWidth', 1);
    legend(hp, lgdtext, 'Location', 'NorthEast');
    legend('boxoff');
    for j = 1:nchannel
        zscore_avg = mean(zscore_bite(:, :, j), 2, 'omitnan');
        zscore_avg_bite(:, i, j) = zscore_avg;
%         [zscore_min, minID] = min(zscore_avg(t >= 0 & t < trange_retrieval(2)));
%         temp = t(t >= 0 & t < trange_retrieval(2));
%         dzscore_bite(i, j) = max(zscore_avg(t > temp(minID)))-zscore_min;
%         
%         peakzscore_bite(i, j) = max(zscore_avg(t > 0));
    end
    
    zscore_chew = result.zscore_chew;
    subplot(1, 3, 3);
    for j = 1:nchannel
        hp(j) = plot_tj_MeanSEM(t, zscore_chew(:, :, j), colors(j, :), colors(j, :), 'Time (s)', 'Z-score', 'aligned to Chew');
    end
    yl = ylim;
    plot([0 0], yl, '--k', 'LineWidth', 1);
    legend(hp, lgdtext, 'Location', 'NorthEast');
    legend('boxoff');
    for j = 1:nchannel
        zscore_avg = mean(zscore_chew(:, :, j), 2, 'omitnan');
        zscore_avg_chew(:, i, j) = zscore_avg;
%         [zscore_min, minID] = min(zscore_avg(t >= 0 & t < trange_retrieval(2)));
%         temp = t(t >= 0 & t < trange_retrieval(2));
%         dzscore_bite(i, j) = max(zscore_avg(t > temp(minID)))-zscore_min;
%         
%         peakzscore_bite(i, j) = max(zscore_avg(t > 0));
    end
end
figure;
subplot(1, 3, 1);
for j = 1:nchannel
    hp(j) = plot_tj_MeanSEM(t, zscore_avg_handadj(:, :, j), colors(j, :), colors(j, :), 'Time (s)', 'Z-score', 'aligned to Hand Adjustment');
end
yl = ylim;
plot([0 0], yl, '--k', 'LineWidth', 1);
legend(hp, lgdtext, 'Location', 'NorthEast');
legend('boxoff');

subplot(1, 3, 2);
for j = 1:nchannel
    hp(j) = plot_tj_MeanSEM(t, zscore_avg_bite(:, :, j), colors(j, :), colors(j, :), 'Time (s)', 'Z-score', 'aligned to Bite');
end
yl = ylim;
plot([0 0], yl, '--k', 'LineWidth', 1);
legend(hp, lgdtext, 'Location', 'NorthEast');
legend('boxoff');

subplot(1, 3, 3);
for j = 1:nchannel
    hp(j) = plot_tj_MeanSEM(t, zscore_avg_chew(:, :, j), colors(j, :), colors(j, :), 'Time (s)', 'Z-score', 'aligned to Chew');
end
yl = ylim;
plot([0 0], yl, '--k', 'LineWidth', 1);
legend(hp, lgdtext, 'Location', 'NorthEast');
legend('boxoff');

cd(curpwd);

%%
n = 6;
figure;
disp('dZ-score hand adjustment');
paired_plot(result(:, 2), result(:, 1), '\Delta Z-score', [], 'dot', [1 2]);

disp('dZ-score bite');
paired_plot(result(:, 4), result(:, 3), '\Delta Z-score', [], 'dot', [3 4]);

disp('dZ-score chew');
paired_plot(result(:, 6), result(:, 5), '\Delta Z-score', [], 'dot', [5 6]);

xlim([0.5 n+0.5]);
set(gca, 'XTick', 1:n, 'XTickLabel', {'aCFA', 'RFO', 'aCFA', 'RFO', 'aCFA', 'RFO'});

%%
% for the comparison with 1 mm pasta lick behavior
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

nchannel = 2;
trange_diff = [-2 2];
colors = [0.3010 0.7450 0.9330; 0.8500 0.3250 0.0980];
hp = zeros(1, nchannel);
for i = 1:nchannel
    lgdtext{i} = ['Channel ' num2str(i)];
end

zscore_retrievalstart_all = cell(1, nchannel);
zscore_lickstart_all = cell(1, nchannel);
zscore_diff = cell(1, nchannel);
SitEnd2RetrievalStartMouth_all = nan(1, N);

for i = 1:N
    temp = load([pathname filename{i}]);
    result = temp.result;
    t = result.t;
%     [~, tzeroID] = min(abs(t));
    tlick = result.tlick;
    zscore_retrievalstart = result.zscore_retrievalstart;
    zscore_lickstart = result.zscore_lickstart;
    SitEnd2RetrievalStartMouth = result.SitEnd2RetrievalStartMouth;
    SitEnd2RetrievalStartMouth_all(i) = mean(SitEnd2RetrievalStartMouth);
    
    for j = 1:nchannel
        figure;
        heatmap4gcamp(tlick(tlick >= -1), zscore_lickstart(tlick >= -1, :, j), ['aligned to Retrieval Start with Mouth: ' lgdtext{j}]);
    end
    
    figure('Name', filename{i});
    for j = 1:nchannel
%         zscore_retrievalstart(:, :, j) = zscore_retrievalstart(:, :, j)-mean(zscore_retrievalstart(tzeroID, :, j), 2);
%         zscore_lickstart(:, :, j) = zscore_lickstart(:, :, j)-mean(zscore_lickstart(tzeroID, :, j), 2);
        
        subplot(1, nchannel, j);
        hp(1) = plot_tj_MeanSEM(tlick, zscore_lickstart(:, :, j), colors(1, :), colors(1, :), 'Time (s)', 'Z-score', 'aligned to Retrieval Start with Mouth');
        hp(2) = plot_tj_MeanSEM(t, zscore_retrievalstart(:, :, j), colors(2, :), colors(2, :), 'Time (s)', 'Z-score', 'aligned to Retrieval Start with Mouth');
        yl = ylim;
        plot([0 0], yl, '--k', 'LineWidth', 1);
        plot(ones(1, 2)*SitEnd2RetrievalStartMouth_all(i), yl, '--k', 'LineWidth', 1);
        legend(hp, {'1 mm pasta', '15 mm pasta'}, 'Location', 'NorthEast');
        legend('boxoff');
        
        if i == 1
            zscore_retrievalstart_all{j} = [];
            zscore_lickstart_all{j} = [];
            zscore_diff{j} = [];
        end
        mean_retrievalstart = mean(zscore_retrievalstart(:, :, j), 2, 'omitnan');
        mean_retrievalstart = mean_retrievalstart(t >= trange_diff(1) & t <= trange_diff(2));
        mean_lickstart = mean(zscore_lickstart(:, :, j), 2, 'omitnan');
        mean_lickstart = mean_lickstart(tlick >= trange_diff(1) & tlick <= trange_diff(2));
        zscore_retrievalstart_all{j} = [zscore_retrievalstart_all{j} mean_retrievalstart];
        zscore_lickstart_all{j} = [zscore_lickstart_all{j} mean_lickstart];
        zscore_diff{j} = [zscore_diff{j} mean_retrievalstart-mean_lickstart];
    end
end
figure;
subplot(1, 3, 1);
hp(1) = plot_tj_MeanSEM(t(t >= trange_diff(1) & t <= trange_diff(2)), zscore_lickstart_all{1}, colors(1, :), colors(1, :), 'Time (s)', 'Z-score', 'aligned to Retrieval Start with Mouth');
hp(2) = plot_tj_MeanSEM(t(t >= trange_diff(1) & t <= trange_diff(2)), zscore_retrievalstart_all{1}, colors(2, :), colors(2, :), 'Time (s)', 'Z-score', 'aligned to Retrieval Start with Mouth');
yl = ylim;
plot([0 0], yl, '--k', 'LineWidth', 1);
plot(ones(1, 2)*mean(SitEnd2RetrievalStartMouth_all), yl, '--k', 'LineWidth', 1);
set(gca, 'ButtonDownFcn', @extract_figure);
legend(hp, {'1 mm pasta', '15 mm pasta'}, 'Location', 'NorthEast');
legend('boxoff');

subplot(1, 3, 2);
hp(1) = plot_tj_MeanSEM(t(t >= trange_diff(1) & t <= trange_diff(2)), zscore_lickstart_all{2}, colors(1, :), colors(1, :), 'Time (s)', 'Z-score', 'aligned to Retrieval Start with Mouth');
hp(2) = plot_tj_MeanSEM(t(t >= trange_diff(1) & t <= trange_diff(2)), zscore_retrievalstart_all{2}, colors(2, :), colors(2, :), 'Time (s)', 'Z-score', 'aligned to Retrieval Start with Mouth');
yl = ylim;
plot([0 0], yl, '--k', 'LineWidth', 1);
plot(ones(1, 2)*mean(SitEnd2RetrievalStartMouth_all), yl, '--k', 'LineWidth', 1);
set(gca, 'ButtonDownFcn', @extract_figure);
legend(hp, {'1 mm pasta', '15 mm pasta'}, 'Location', 'NorthEast');
legend('boxoff');

subplot(1, 3, 3);
for i = 1:nchannel
    hp(i) = plot_tj_MeanSEM(t(t >= trange_diff(1) & t <= trange_diff(2)), zscore_diff{i}, colors(i, :), colors(i, :), 'Time (s)', 'dZ-score', 'aligned to Retrieval Start with Mouth');
end
yl = ylim;
plot([0 0], yl, '--k', 'LineWidth', 1);
plot(ones(1, 2)*mean(SitEnd2RetrievalStartMouth_all), yl, '--k', 'LineWidth', 1);
set(gca, 'ButtonDownFcn', @extract_figure);
legend(hp, lgdtext, 'Location', 'NorthEast');
legend('boxoff');

cd(curpwd);

%%
% for the comparison with pasta bite device
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

nchannel = 2;
trange = [-2 5];
colors = [0.3010 0.7450 0.9330; 0.8500 0.3250 0.0980];
hp = zeros(1, nchannel);
for i = 1:nchannel
    lgdtext{i} = ['Channel ' num2str(i)];
end

zscore_firstbite_all = cell(1, nchannel);

for i = 1:N
    temp = load([pathname filename{i}]);
    result = temp.result;
    tbite = result.tbite;
    zscore_firstbite = result.zscore_firstbite;
    
    for j = 1:nchannel
        figure;
        heatmap4gcamp(tbite(tbite >= trange(1) & tbite <= trange(2)), zscore_firstbite(tbite >= trange(1) & tbite <= trange(2), :, j), ['aligned to first bite: ' lgdtext{j}]);
    end
    
    figure('Name', filename{i});
    for j = 1:nchannel
        hp(j) = plot_tj_MeanSEM(tbite(tbite >= trange(1) & tbite <= trange(2)), zscore_firstbite(tbite >= trange(1) & tbite <= trange(2), :, j),...
            colors(j, :), colors(j, :), 'Time (s)', 'Z-score', 'aligned to first bite');
        
        if i == 1
            zscore_firstbite_all{j} = [];
        end
        mean_firstbite = mean(zscore_firstbite(:, :, j), 2, 'omitnan');
        mean_firstbite = mean_firstbite(tbite >= trange(1) & tbite <= trange(2));
        zscore_firstbite_all{j} = [zscore_firstbite_all{j} mean_firstbite];
    end
    yl = ylim;
    plot([0 0], yl, '--k', 'LineWidth', 1);
    legend(hp, lgdtext, 'Location', 'NorthEast');
    legend('boxoff');
end
figure;
for i = 1:nchannel
    hp(i) = plot_tj_MeanSEM(tbite(tbite >= trange(1) & tbite <= trange(2)), zscore_firstbite_all{i}, colors(i, :), colors(i, :), 'Time (s)', 'Z-score', 'aligned to first bite');
end
yl = ylim;
plot([0 0], yl, '--k', 'LineWidth', 1);
legend(hp, lgdtext, 'Location', 'NorthEast');
legend('boxoff');

cd(curpwd);

%%
% for step crossing across mice
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

nchannel = 2;
hp = zeros(1, nchannel);
for i = 1:nchannel
    lgdtext{i} = ['Channel ' num2str(i)];
end

zscore_avg_stepcross = [];

dzscore_stepcross = nan(N, nchannel);
peakzscore_stepcross = nan(N, nchannel);

retrievalstart2stepcross_all = nan(1, N);

for i = 1:N
    temp = load([pathname filename{i}]);
    result = temp.result;
    t = result.t;
    zscore_stepcross = result.zscore_stepcross;
    retrievalstart2stepcross_all(i) = mean(result.retrievalstart2stepcross);
    figure('Name', filename{i});
    for j = 1:nchannel
        hp(j) = plot_tj_MeanSEM(t, zscore_stepcross(:, :, j), [1-1/nchannel*j 1-1/nchannel*j 1-1/nchannel*j], [1-1/nchannel*j 1-1/nchannel*j 1-1/nchannel*j], 'Time (s)', 'Z-score', 'Aligned to step cross');
    end
    yl = ylim;
    plot([0 0], yl, '--k', 'LineWidth', 1);
    legend(hp, lgdtext, 'Location', 'NorthEast');
    legend('boxoff');
    for j = 1:nchannel
        zscore_avg = mean(zscore_stepcross(:, :, j), 2, 'omitnan');
        zscore_avg_stepcross(:, i, j) = zscore_avg;
    end
end
figure;
for j = 1:nchannel
    hp(j) = plot_tj_MeanSEM(t, zscore_avg_stepcross(:, :, j), [1-1/nchannel*j 1-1/nchannel*j 1-1/nchannel*j], [1-1/nchannel*j 1-1/nchannel*j 1-1/nchannel*j], 'Time (s)', 'Z-score', 'Aligned to step cross');
end
yl = ylim;
plot([0 0], yl, '--k', 'LineWidth', 1);
plot(ones(1, 2)*mean(retrievalstart2stepcross_all), yl, '--k', 'LineWidth', 1);
legend(hp, lgdtext, 'Location', 'NorthEast');
legend('boxoff');

cd(curpwd);