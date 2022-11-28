%%
N = 2;
result.cc_all = [];
for i = 1:N
    eval(['result.cc_all = [result.cc_all result' num2str(i) '.cc_all];']);
end

%% GCaMP signal and hand-to-nose distance correlation coefficient
clc;
ccstep = 1; % number of data point, default 1, shift 1 data point each time
shifttime = 120; % half of shift time, default 120, shift 1 sec for both left and right sides
FrameRate = 120;
nchannel = 2;
lgdtext = {'Channel 1', 'Channel 2'};
t_cc = (-shifttime:1:shifttime)'*ccstep/FrameRate;

curpwd = pwd;
try
   cd(pathname); 
end
[filename, pathname] = uigetfile('*.mat', 'Pick all of the PlexinD1 data set', 'MultiSelect', 'on');
if isequal(filename, 0)
    cd(curpwd);
    return;
end

N = numel(filename);

cc_peakP = nan(1, N);
t_peakP = nan(1, N);

for i = 1:N
    temp = load([pathname filename{i}]);
    result = temp.result;
    cc_avg = mean(result.cc_all(:, :, 2), 2);
    [cc_peakP(i), maxID] = max(cc_avg);
    t_peakP(i) = t_cc(maxID);
    
    figure;
    subplot(1, 2, 1);
    plot_tj_individuals(t_cc, result.cc_all(:, :, 2), [0.75 0.75 0.75], [0 0 0], 'Lag (s)', 'Correlation coefficient', ['Z-score with hand-to-nose distance: ' lgdtext{2}]);
    
    subplot(1, 2, 2);
    hp = plot_tj_MeanSEM(t_cc, result.cc_all(:, :, 2), [1 0.41 0.16], [0.8 0 0], 'Lag (s)', 'Correlation coefficient', lgdtext{2});
    yl = ylim;
    line([0 0], yl, 'Color', [0 0 0], 'LineStyle', '--', 'LineWidth', 1);
    legend(hp, 'Z-score with hand-to-nose distance');
    legend('boxoff');
end

try
   cd(pathname); 
end
[filename, pathname] = uigetfile('*.mat', 'Pick all of the Fezf2 data set', 'MultiSelect', 'on');
if isequal(filename, 0)
    cd(curpwd);
    return;
end

N = numel(filename);

cc_peakF = nan(1, N);
t_peakF = nan(1, N);

for i = 1:N
    temp = load([pathname filename{i}]);
    result = temp.result;
    cc_avg = mean(result.cc_all(:, :, 2), 2);
    [cc_peakF(i), maxID] = max(cc_avg);
    t_peakF(i) = t_cc(maxID);
    
    figure;
    subplot(1, 2, 1);
    plot_tj_individuals(t_cc, result.cc_all(:, :, 2), [0.75 0.75 0.75], [0 0 0], 'Lag (s)', 'Correlation coefficient', ['Z-score with hand-to-nose distance: ' lgdtext{2}]);
    
    subplot(1, 2, 2);
    hp = plot_tj_MeanSEM(t_cc, result.cc_all(:, :, 2), [1 0.41 0.16], [0.8 0 0], 'Lag (s)', 'Correlation coefficient', lgdtext{2});
    yl = ylim;
    line([0 0], yl, 'Color', [0 0 0], 'LineStyle', '--', 'LineWidth', 1);
    legend(hp, 'Z-score with hand-to-nose distance');
    legend('boxoff');
end

violin_boxplot(cc_peakF, cc_peakP, 'Peak correlation coefficient');
violin_boxplot(t_peakF, t_peakP, 'Peak time');

cd(curpwd);