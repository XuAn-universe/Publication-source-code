%%
% preprocess raw data
AnimalFolder = 'G:\Fiber Photometry for CellReadr\20220405_Ctip2_Gcamp2';
ExpFolder = [];
curpwd = pwd;
cd(AnimalFolder);
DirList = dir;
totalfile = numel(DirList);
n = 0;
workbar(0, 'Computing Ongoing...', 'Progress');
for i = 3:totalfile
    if numel(DirList(i).name) < 3
        continue;
    end
    if strcmp(DirList(i).name(1:3), 'exp')
        n = n+1;
        ExpFolder{n} = [AnimalFolder '\' DirList(i).name];
        Raw_FPData_Preprocessing([ExpFolder{n} '\signal0.csv'],...
            [ExpFolder{n} '\trigger0.csv'],...
            ExpFolder{n}, [], 1, 1e12, 0, 5, 0, 1e10);
        close all;
    end
    workbar(i/totalfile, [num2str(i) '/' num2str(totalfile)], 'Progress');
end
pause(0.5);
cd(curpwd);

%%
% align signal to stimulus onset for single experiment 
exp = 9;
pre = 1;
post = 4;
channel = 2;
result = load([ExpFolder{exp} '\FPData.mat']);
result = result.zsignal_all;
ntrial = size(result, 2);
signal_all = [];
signal_random_all = [];
figure;
for i = 1:ntrial
    try
        subplot(4, 5, i);
        t = result{channel, i}(:, 1);
        signal = result{channel, i}(:, 2);
        plot(t, signal, '-k');
        hold on;
        plot([t(1) t(end)], [0 0], '--k');
        yl = ylim;
        plot([0 0], yl, '--r', 'LineWidth', 1);
        plot([1 1], yl, '--r', 'LineWidth', 1);
        box off;
        set(gca, 'ButtonDownFcn', @extract_figure);
    end
    
    [t_aligned, zscore_signal, zscore_signal_random] = AlignSignal2Event(t, signal, 0, pre, post);
    try
        signal_all = [signal_all zscore_signal];
        signal_random_all = [signal_random_all zscore_signal_random];
    catch
        if numel(zscore_signal) < size(signal_all, 1)
            zscore_signal(end+1) = zscore_signal(end);
            zscore_signal_random(end+1) = zscore_signal_random(end);
            signal_all = [signal_all zscore_signal];
            signal_random_all = [signal_random_all zscore_signal_random];
        end
    end
end

IDexc = [];
signal_all(:, IDexc) = [];
signal_random_all(:, IDexc) = [];

figure;
heatmap4gcamp(t_aligned, signal_all, '');

colors = [0.3010 0.7450 0.9330; 0 0 0];
figure;
plot_tj_MeanSEM(t_aligned, signal_all, colors(1, :), colors(1, :), 'Time (s)', 'Z-score', '');
plot_tj_MeanSEM(t_aligned, signal_random_all, colors(2, :), colors(2, :), 'Time (s)', 'Z-score', '');
yl = ylim;
plot([0 0], yl, '--r', 'LineWidth', 1);
plot([1 1], yl, '--r', 'LineWidth', 1);

%%
% combine signals from different experiment
signal_all = [signal_all8 signal_all9 signal_all4];
signal_random_all = [signal_random_all8 signal_random_all9 signal_random_all4];

figure;
heatmap4gcamp(t_aligned(t_aligned <= 4), signal_all(t_aligned <= 4, :), '');

colors = [0.3010 0.7450 0.9330; 0 0 0];
figure;
plot_tj_MeanSEM(t_aligned(t_aligned <= 4), signal_all(t_aligned <= 4, :), colors(1, :), colors(1, :), 'Time (s)', 'Z-score', '');
plot_tj_MeanSEM(t_aligned(t_aligned <= 4), signal_random_all(t_aligned <= 4, :), colors(2, :), colors(2, :), 'Time (s)', 'Z-score', '');
yl = ylim;
plot([0 0], yl, '--r', 'LineWidth', 1);
plot([1 1], yl, '--r', 'LineWidth', 1);

%%
% get the peak Z-score
signal_avg = mean(signal_all, 2, 'omitnan');
[Zpeak, IDpeak] = max(signal_avg);
figure;
plot(t_aligned, signal_avg, '-k');
hold on;
plot(t_aligned(IDpeak), Zpeak, 'ro')