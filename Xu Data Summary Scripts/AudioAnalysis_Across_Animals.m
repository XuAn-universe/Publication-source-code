%*---------------------------------------------------------------------*
% Huang Lab
% Cold Spring Harbor Laboratory
% Author : Xu An, Sep 2019
% xan@cshl.edu
% Version: 1.0
%*---------------------------------------------------------------------*
%%
[filename, pathname] = uigetfile('*.mat', 'Pick all of the individual data set', 'MultiSelect', 'on');
if isequal(filename, 0)
    return;
end
clc;
N = numel(filename);
resampling = 10000;
firstafter2lastbefore_I = cell(1, N);
firstafter2lastbefore_NI = cell(1, N);
Afirstafter2lastbefore_I = zeros(1, N); % 'A' means average
Afirstafter2lastbefore_NI = zeros(1, N);

totalbites_I = cell(1, N);
totalbites_NI = cell(1, N);
Atotalbites_I = zeros(1, N);
Atotalbites_NI = zeros(1, N);

bite_duration_I = cell(1, N);
bite_duration_NI = cell(1, N);
Abite_duration_I = zeros(1, N);
Abite_duration_NI = zeros(1, N);

bite_intervals_I = cell(1, N);
bite_intervals_NI = cell(1, N);
Abite_intervals_I = zeros(1, N);
Abite_intervals_NI = zeros(1, N);

bite_amplitudes_I = cell(1, N);
bite_amplitudes_NI = cell(1, N);
Abite_amplitudes_I = zeros(1, N);
Abite_amplitudes_NI = zeros(1, N);

firstbite_time_I = cell(1, N);
firstbite_time_NI = cell(1, N);
Afirstbite_time_I = zeros(1, N);
Afirstbite_time_NI = zeros(1, N);
for i = 1:N
    load([pathname filename{i}]);
    if all(isfield(result, {'firstafter2lastbefore_I', 'firstafter2lastbefore_NI'}))
        Afirstafter2lastbefore_I(i) = mean(result.firstafter2lastbefore_I);
        Afirstafter2lastbefore_NI(i) = mean(result.firstafter2lastbefore_NI);
        firstafter2lastbefore_I{i} = result.firstafter2lastbefore_I;
        firstafter2lastbefore_NI{i} = result.firstafter2lastbefore_NI;
    elseif all(isfield(result, {'firstbite_time_I', 'firstbite_time_NI'}))
        Afirstbite_time_I(i) = mean(result.firstbite_time_I);
        Afirstbite_time_NI(i) = mean(result.firstbite_time_NI);
        firstbite_time_I{i} = result.firstbite_time_I;
        firstbite_time_NI{i} = result.firstbite_time_NI;
        if all(isfield(result, {'bite_duration_I', 'bite_duration_NI'}))
            Abite_duration_I(i) = mean(result.bite_duration_I);
            Abite_duration_NI(i) = mean(result.bite_duration_NI);
            bite_duration_I{i} = result.bite_duration_I;
            bite_duration_NI{i} = result.bite_duration_NI;
        end
    end
    Atotalbites_I(i) = mean(result.totalbites_I);
    Atotalbites_NI(i) = mean(result.totalbites_NI);
    totalbites_I{i} = result.totalbites_I;
    totalbites_NI{i} = result.totalbites_NI;
    
    Abite_intervals_I(i) = mean(result.bite_intervals_I);
    Abite_intervals_NI(i) = mean(result.bite_intervals_NI);
    bite_intervals_I{i} = result.bite_intervals_I;
    bite_intervals_NI{i} = result.bite_intervals_NI;
    
    Abite_amplitudes_I(i) = mean(result.bite_amplitudes_I);
    Abite_amplitudes_NI(i) = mean(result.bite_amplitudes_NI);
    bite_amplitudes_I{i} = result.bite_amplitudes_I;
    bite_amplitudes_NI{i} = result.bite_amplitudes_NI;
end

% permutation test
figure;
subplot(2, 2, 1);
if all(isfield(result, {'firstafter2lastbefore_I', 'firstafter2lastbefore_NI'}))
    permutation_test(firstafter2lastbefore_NI, firstafter2lastbefore_I, 'both', resampling, 'difference of firstafter2lastbefore (s)');
elseif all(isfield(result, {'firstbite_time_I', 'firstbite_time_NI'}))
    permutation_test(firstbite_time_NI, firstbite_time_I, 'both', resampling, 'difference of firstbite time (s)');
end
subplot(2, 2, 2);
permutation_test(totalbites_NI, totalbites_I, 'both', resampling, 'difference of totalbites');
subplot(2, 2, 3);
permutation_test(bite_intervals_NI, bite_intervals_I, 'both', resampling, 'difference of bite interval (s)');
subplot(2, 2, 4);
permutation_test(bite_amplitudes_NI, bite_amplitudes_I, 'both', resampling, 'difference of bite amplitude');


if all(isfield(result, {'bite_duration_I', 'bite_duration_NI'}))
    figure;
    permutation_test(bite_duration_NI, bite_duration_I, 'both', resampling, 'difference of total duration of all bites (s)');
end

% bootstrapping test
figure;
subplot(2, 2, 1);
if all(isfield(result, {'firstafter2lastbefore_I', 'firstafter2lastbefore_NI'}))
    bootstrapping_test(firstafter2lastbefore_NI, firstafter2lastbefore_I, 'both', resampling, 'difference of firstafter2lastbefore (s)');
elseif all(isfield(result, {'firstbite_time_I', 'firstbite_time_NI'}))
    bootstrapping_test(firstbite_time_NI, firstbite_time_I, 'both', resampling, 'difference of firstbite time (s)');
end
subplot(2, 2, 2);
bootstrapping_test(totalbites_NI, totalbites_I, 'both', resampling, 'difference of totalbites');
subplot(2, 2, 3);
bootstrapping_test(bite_intervals_NI, bite_intervals_I, 'both', resampling, 'difference of bite interval (s)');
subplot(2, 2, 4);
bootstrapping_test(bite_amplitudes_NI, bite_amplitudes_I, 'both', resampling, 'difference of bite amplitude');

if all(isfield(result, {'bite_duration_I', 'bite_duration_NI'}))
    figure;
    bootstrapping_test(bite_duration_NI, bite_duration_I, 'both', resampling, 'difference of total duration of all bites (s)');
end

if all(isfield(result, {'firstafter2lastbefore_I', 'firstafter2lastbefore_NI'}))
    firstafter2lastbefore_NI = cell2vector(firstafter2lastbefore_NI);
    firstafter2lastbefore_I = cell2vector(firstafter2lastbefore_I);
elseif all(isfield(result, {'firstbite_time_I', 'firstbite_time_NI'}))
    firstbite_time_NI = cell2vector(firstbite_time_NI);
    firstbite_time_I = cell2vector(firstbite_time_I);
    if all(isfield(result, {'bite_duration_I', 'bite_duration_NI'}))
        bite_duration_NI = cell2vector(bite_duration_NI);
        bite_duration_I = cell2vector(bite_duration_I);
    end
end
totalbites_NI = cell2vector(totalbites_NI);
totalbites_I = cell2vector(totalbites_I);
bite_intervals_NI = cell2vector(bite_intervals_NI);
bite_intervals_I = cell2vector(bite_intervals_I);
bite_amplitudes_NI = cell2vector(bite_amplitudes_NI);
bite_amplitudes_I = cell2vector(bite_amplitudes_I);

% violin plots
if all(isfield(result, {'firstafter2lastbefore_I', 'firstafter2lastbefore_NI'}))
    violin_boxplot(firstafter2lastbefore_NI, firstafter2lastbefore_I, 'Interval (s)');
elseif all(isfield(result, {'firstbite_time_I', 'firstbite_time_NI'}))
    violin_boxplot(firstbite_time_NI, firstbite_time_I, 'Interval (s)');
end
violin_boxplot(totalbites_NI, totalbites_I, 'Bites');
violin_boxplot(bite_intervals_NI, bite_intervals_I, 'Bite Interval (s)');
violin_boxplot(bite_amplitudes_NI, bite_amplitudes_I, 'Normalized Sound Amplitude');

if all(isfield(result, {'bite_duration_I', 'bite_duration_NI'}))
    violin_boxplot(bite_duration_NI, bite_duration_I, 'total duration of all bites (s)');
end

% bar plots
linecolor = [0.8 0.8 0.8];
figure;
subplot(1, 4, 1);
hold on;
if all(isfield(result, {'firstafter2lastbefore_I', 'firstafter2lastbefore_NI'}))
    for i = 1:N
        plot([1 2], [Afirstafter2lastbefore_NI(i) Afirstafter2lastbefore_I(i)], '-o', 'Color', linecolor, 'MarkerFaceColor', linecolor, 'MarkerEdgeColor', 'none', 'ButtonDownFcn', @displayfilename, 'UserData', filename{i});
    end
    plot([1 2], [mean(Afirstafter2lastbefore_NI) mean(Afirstafter2lastbefore_I)], '-ok', 'MarkerFaceColor', [0 0 0], 'MarkerEdgeColor', 'none', 'LineWidth', 2);
    set(gca, 'XTick', [1 2], 'XTickLabel', {'Control', 'Inhibition'}, 'XLim', [0.4 2.6], 'FontSize', 12, 'ButtonDownFcn', @extract_figure);
    ylabel('Interval (s)');
    yl = ylim;
    ylim([0 yl(2)]);
    fprintf('\nPaired t test:\n');
    [~, p] = ttest(Afirstafter2lastbefore_NI, Afirstafter2lastbefore_I, 'tail', 'both');
    disp(['firstafter2lastbefore (both): p = ' num2str(p)]);
    [H, ~, ~] = swtest(Afirstafter2lastbefore_NI-Afirstafter2lastbefore_I);
    if H
        disp('Data didn''t pass normality test');
    else
        disp('Data passed normality test');
    end
    fprintf('\nSignrank test:\n');
    p = signrank(Afirstafter2lastbefore_NI, Afirstafter2lastbefore_I, 'tail', 'both');
    disp(['firstafter2lastbefore (both): p = ' num2str(p)]);
elseif all(isfield(result, {'firstbite_time_I', 'firstbite_time_NI'}))
    for i = 1:N
        plot([1 2], [Afirstbite_time_NI(i) Afirstbite_time_I(i)], '-o', 'Color', linecolor, 'MarkerFaceColor', linecolor, 'MarkerEdgeColor', 'none', 'ButtonDownFcn', @displayfilename, 'UserData', filename{i});
    end
    plot([1 2], [mean(Afirstbite_time_NI) mean(Afirstbite_time_I)], '-ok', 'MarkerFaceColor', [0 0 0], 'MarkerEdgeColor', 'none', 'LineWidth', 2);
    set(gca, 'XTick', [1 2], 'XTickLabel', {'Control', 'Inhibition'}, 'XLim', [0.4 2.6], 'FontSize', 12, 'ButtonDownFcn', @extract_figure);
    ylabel('Interval (s)');
    yl = ylim;
    ylim([0 yl(2)]);
    fprintf('\nPaired t test:\n');
    [~, p] = ttest(Afirstbite_time_NI, Afirstbite_time_I, 'tail', 'both');
    disp(['First Bite (both): p = ' num2str(p)]);
    [H, ~, ~] = swtest(Afirstbite_time_NI-Afirstbite_time_I);
    if H
        disp('Data didn''t pass normality test');
    else
        disp('Data passed normality test');
    end
    fprintf('\nSignrank test:\n');
    p = signrank(Afirstbite_time_NI, Afirstbite_time_I, 'tail', 'both');
    disp(['First Bite (both): p = ' num2str(p)]);
end

subplot(1, 4, 2);
hold on;
for i = 1:N
    plot([1 2], [Atotalbites_NI(i) Atotalbites_I(i)], '-o', 'Color', linecolor, 'MarkerFaceColor', linecolor, 'MarkerEdgeColor', 'none', 'ButtonDownFcn', @displayfilename, 'UserData', filename{i});
end
plot([1 2], [mean(Atotalbites_NI) mean(Atotalbites_I)], '-ok', 'MarkerFaceColor', [0 0 0], 'MarkerEdgeColor', 'none', 'LineWidth', 2);
set(gca, 'XTick', [1 2], 'XTickLabel', {'Control', 'Inhibition'}, 'XLim', [0.4 2.6], 'FontSize', 12, 'ButtonDownFcn', @extract_figure);
ylabel('Bites');
yl = ylim;
ylim([0 yl(2)]);
fprintf('\nPaired t test:\n');
[~, p] = ttest(Atotalbites_NI, Atotalbites_I, 'tail', 'both');
disp(['totalbites (both): p = ' num2str(p)]);
[H, ~, ~] = swtest(Atotalbites_NI-Atotalbites_I);
if H
    disp('Data didn''t pass normality test');
else
    disp('Data passed normality test');
end
fprintf('\nSignrank test:\n');
p = signrank(Atotalbites_NI, Atotalbites_I, 'tail', 'both');
disp(['totalbites (both): p = ' num2str(p)]);

try
    subplot(1, 4, 3);
    hold on;
    for i = 1:N
        plot([1 2], [Abite_intervals_NI(i) Abite_intervals_I(i)], '-o', 'Color', linecolor, 'MarkerFaceColor', linecolor, 'MarkerEdgeColor', 'none', 'ButtonDownFcn', @displayfilename, 'UserData', filename{i});
    end
    plot([1 2], [mean(Abite_intervals_NI) mean(Abite_intervals_I)], '-ok', 'MarkerFaceColor', [0 0 0], 'MarkerEdgeColor', 'none', 'LineWidth', 2);
    set(gca, 'XTick', [1 2], 'XTickLabel', {'Control', 'Inhibition'}, 'XLim', [0.4 2.6], 'FontSize', 12, 'ButtonDownFcn', @extract_figure);
    ylabel('Bite Interval (s)');
    yl = ylim;
    ylim([0 yl(2)]);
    fprintf('\nPaired t test:\n');
    [~, p] = ttest(Abite_intervals_NI, Abite_intervals_I, 'tail', 'both');
    disp(['bite interval (both): p = ' num2str(p)]);
    [H, ~, ~] = swtest(Abite_intervals_NI-Abite_intervals_I);
    if H
        disp('Data didn''t pass normality test');
    else
        disp('Data passed normality test');
    end
    fprintf('\nSignrank test:\n');
    p = signrank(Abite_intervals_NI, Abite_intervals_I, 'tail', 'both');
    disp(['bite interval (both): p = ' num2str(p)]);
end

try
    subplot(1, 4, 4);
    hold on;
    for i = 1:N
        plot([1 2], [Abite_amplitudes_NI(i) Abite_amplitudes_I(i)], '-o', 'Color', linecolor, 'MarkerFaceColor', linecolor, 'MarkerEdgeColor', 'none', 'ButtonDownFcn', @displayfilename, 'UserData', filename{i});
    end
    plot([1 2], [mean(Abite_amplitudes_NI) mean(Abite_amplitudes_I)], '-ok', 'MarkerFaceColor', [0 0 0], 'MarkerEdgeColor', 'none', 'LineWidth', 2);
    set(gca, 'XTick', [1 2], 'XTickLabel', {'Control', 'Inhibition'}, 'XLim', [0.4 2.6], 'FontSize', 12, 'ButtonDownFcn', @extract_figure);
    ylabel('Normalized Sound Amplitude');
    ylim([0 1]);
    fprintf('\nPaired t test:\n');
    [~, p] = ttest(Abite_amplitudes_NI, Abite_amplitudes_I, 'tail', 'both');
    disp(['bite amplitude (both): p = ' num2str(p)]);
    [H, ~, ~] = swtest(Abite_amplitudes_NI-Abite_amplitudes_I);
    if H
        disp('Data didn''t pass normality test');
    else
        disp('Data passed normality test');
    end
    fprintf('\nSignrank test:\n');
    p = signrank(Abite_amplitudes_NI, Abite_amplitudes_I, 'tail', 'both');
    disp(['bite amplitude (both): p = ' num2str(p)]);
end

if all(isfield(result, {'bite_duration_I', 'bite_duration_NI'}))
    figure;
    hold on;
    for i = 1:N
        plot([1 2], [Abite_duration_NI(i) Abite_duration_I(i)], '-o', 'Color', linecolor, 'MarkerFaceColor', linecolor, 'MarkerEdgeColor', 'none', 'ButtonDownFcn', @displayfilename, 'UserData', filename{i});
    end
    plot([1 2], [mean(Abite_duration_NI) mean(Abite_duration_I)], '-ok', 'MarkerFaceColor', [0 0 0], 'MarkerEdgeColor', 'none', 'LineWidth', 2);
    set(gca, 'XTick', [1 2], 'XTickLabel', {'Control', 'Inhibition'}, 'XLim', [0.4 2.6], 'FontSize', 12, 'ButtonDownFcn', @extract_figure);
    ylabel('total duration of all bites (s)');
    yl = ylim;
    ylim([0 yl(2)]);
    fprintf('\nPaired t test:\n');
    [~, p] = ttest(Abite_duration_NI, Abite_duration_I, 'tail', 'both');
    disp(['total duration of all bites (both): p = ' num2str(p)]);
    [H, ~, ~] = swtest(Abite_duration_NI-Abite_duration_I);
    if H
        disp('Data didn''t pass normality test');
    else
        disp('Data passed normality test');
    end
    fprintf('\nSignrank test:\n');
    p = signrank(Abite_duration_NI, Abite_duration_I, 'tail', 'both');
    disp(['total duration of all bites (both): p = ' num2str(p)]);
end

% histograms & cdf
try
    [~, edges] = histcounts([bite_intervals_I bite_intervals_NI]);
    figure;
    subplot(2, 2, 1);
    histogram(bite_intervals_NI, edges, 'FaceColor', [0 0 0], 'FaceAlpha', 1, 'Normalization', 'probability');
    hold on;
    histogram(bite_intervals_I, edges, 'EdgeColor', [0 1 0], 'FaceColor', 'none', 'LineWidth', 1, 'Normalization', 'probability');
    set(gca, 'TickLength', [0 0], 'FontSize', 12);
    box off;
    xlabel('Bite Interval (s)');
    ylabel('Probability');
    subplot(2, 2, 3);
    hold on;
    h = cdfplot(bite_intervals_NI);
    set(h, 'Color', [0 0 0], 'LineWidth', 1.5);
    if ~isempty(bite_intervals_I)
        h = cdfplot(bite_intervals_I);
        set(h, 'Color', [0 1 0], 'LineWidth', 1.5);
    end
    set(gca, 'GridLineStyle', ':', 'FontSize', 12);
    xlabel('Bite Interval (s)');
    ylabel('Cumulative Probability');
    fprintf('\n');
    if ~isempty(bite_intervals_I)
        [~, p] = kstest2(bite_intervals_NI, bite_intervals_I);
        disp(['bite_intervals (KStest): p = ' num2str(p)]);
    end
end

% [~, edges] = histcounts([bite_amplitudes_I bite_amplitudes_NI]);
try
    edges = 0:0.1:1;
    subplot(2, 2, 2);
    histogram(bite_amplitudes_NI, edges, 'FaceColor', [0 0 0], 'FaceAlpha', 1, 'Normalization', 'probability');
    hold on;
    histogram(bite_amplitudes_I, edges, 'EdgeColor', [0 1 0], 'FaceColor', 'none', 'LineWidth', 1, 'Normalization', 'probability');
    set(gca, 'TickLength', [0 0], 'FontSize', 12);
    box off;
    xlabel('Normalized Sound Amplitude');
    ylabel('Probability');
    subplot(2, 2, 4);
    hold on;
    h = cdfplot(bite_amplitudes_NI);
    set(h, 'Color', [0 0 0], 'LineWidth', 1.5);
    if ~isempty(bite_amplitudes_I)
        h = cdfplot(bite_amplitudes_I);
        set(h, 'Color', [0 1 0], 'LineWidth', 1.5);
    end
    set(gca, 'GridLineStyle', ':', 'FontSize', 12);
    xlabel('Normalized Sound Amplitude');
    ylabel('Cumulative Probability');
    fprintf('\n');
    if ~isempty(bite_amplitudes_I)
        [~, p] = kstest2(bite_amplitudes_NI, bite_amplitudes_I);
        disp(['bite_amplitudes (KStest): p = ' num2str(p)]);
    end
end

%%
% combine result from different sessions
names = fieldnames(result1);
for i = 1:numel(names)
    eval(['result.' names{i} ' = [result1.' names{i} ' result2.' names{i} '];']);
end

%%
[filename, pathname] = uigetfile('*.mat', 'Pick all of the individual data set', 'MultiSelect', 'on');
if isequal(filename, 0)
    return;
end
FrameRate = 120;
N = numel(filename);
totalbite = nan(1, N);
bite_interval_all = [];
for i = 1:N
    load([pathname filename{i}]);
    totalbite(i) = mean(result.totalbites);
    bite_interval_all = [bite_interval_all result.bite_intervals];
end
figure;
singlebar_plot(totalbite, 1, '', 'Bite');

figure;
plothistogram(bite_interval_all, 0:1/FrameRate*2:max(bite_interval_all), [0 0 0], [0 0 0], 'Bite interval (s)')

function vectordata = cell2vector(celldata)
vectordata = [];
for i = 1:numel(celldata)
    vectordata = [vectordata celldata{i}];
end
end

function displayfilename(src, eventdata)
htext = text(eventdata.IntersectionPoint(1), eventdata.IntersectionPoint(2), src.UserData(1:end-4));
pause(2);
try
    delete(htext)
    clear htext;
end
end