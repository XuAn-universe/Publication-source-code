%*---------------------------------------------------------------------*
% Huang Lab
% Duke University
% Author : Xu An, Aug 2022
% xu.an@duke.edu
% Version: 1.0
%*---------------------------------------------------------------------*
%%
N = 2;
result.dzwithdraw_NI = [];
result.dzspeedwithdraw_NI = [];
result.dzwithdraw_I = [];
result.dzspeedwithdraw_I = [];
for i = 1:N
    eval(['result.dzwithdraw_NI = [result.dzwithdraw_NI result' num2str(i) '.dzwithdraw_NI{1}];']);
    eval(['result.dzspeedwithdraw_NI = [result.dzspeedwithdraw_NI result' num2str(i) '.dzspeedwithdraw_NI{1}];']);
    eval(['result.dzwithdraw_I = [result.dzwithdraw_I result' num2str(i) '.dzwithdraw_I{1}];']);
    eval(['result.dzspeedwithdraw_I = [result.dzspeedwithdraw_I result' num2str(i) '.dzspeedwithdraw_I{1}];']);
end

%% hand-to-nose distance aligned to withdraw
curpwd = pwd;
try
   cd(pathname); 
end
[filename, pathname] = uigetfile('*.mat', 'Pick all of the data set', 'MultiSelect', 'on');
if isequal(filename, 0)
    cd(curpwd);
    return;
end

FrameRate = 120;
trange = 1;
trange = trange*FrameRate;
t = (-trange:1:trange)'/FrameRate;

N = numel(filename);
dzwithdraw_NI = nan(2*trange+1, N);
dzwithdraw_I = nan(2*trange+1, N);
dzspeedwithdraw_NI = nan(2*trange+1, N);
dzspeedwithdraw_I = nan(2*trange+1, N);
for i = 1:N
    temp = load([pathname filename{i}]);
    result = temp.result;
    if iscell(result.dzwithdraw_NI)
        dzwithdraw_NI(:, i) = mean(result.dzwithdraw_NI{1}, 2, 'omitnan');
    else
        dzwithdraw_NI(:, i) = mean(result.dzwithdraw_NI, 2, 'omitnan');
    end
    if iscell(result.dzwithdraw_I)
        dzwithdraw_I(:, i) = mean(result.dzwithdraw_I{1}, 2, 'omitnan');
    else
        dzwithdraw_I(:, i) = mean(result.dzwithdraw_I, 2, 'omitnan');
    end
    
    if iscell(result.dzspeedwithdraw_NI)
        dzspeedwithdraw_NI(:, i) = mean(result.dzspeedwithdraw_NI{1}, 2, 'omitnan');
    else
        dzspeedwithdraw_NI(:, i) = mean(result.dzspeedwithdraw_NI, 2, 'omitnan');
    end
    if iscell(result.dzspeedwithdraw_I)
        dzspeedwithdraw_I(:, i) = mean(result.dzspeedwithdraw_I{1}, 2, 'omitnan');
    else
        dzspeedwithdraw_I(:, i) = mean(result.dzspeedwithdraw_I, 2, 'omitnan');
    end
end
dzwithdraw_diff = dzwithdraw_I-dzwithdraw_NI;
figure;
plot_tj_individuals(t, dzwithdraw_NI, [0.8 0.8 0.8], [0 0 0], 'Time (s)', 'dZ (mm)', '', filename);
plot_tj_individuals(t, dzwithdraw_I, [0 1 1], [0 0 1], 'Time (s)', 'dZ (mm)', '', filename);
yl = ylim;
plot([0 0], yl, '--k');
figure;
plot_tj_acrosstrial(t, dzwithdraw_NI, [0 0 0], [0 0 0], 'dZ (mm)');
plot_tj_acrosstrial(t, dzwithdraw_I, [0 0 1], [0 0 1], 'dZ (mm)');
yl = ylim;
plot([0 0], yl, '--k');
figure;
plot_tj_individuals(t, dzwithdraw_diff, [0.8 0.8 0.8], [0 0 0], 'Time (s)', 'ddZ (mm)', '', filename);
yl = ylim;
plot([0 0], yl, '--k');
plot([t(1) t(end)], [0 0], '--k');

dzspeedwithdraw_diff = dzspeedwithdraw_I-dzspeedwithdraw_NI;
figure;
plot_tj_individuals(t, dzspeedwithdraw_NI, [0.8 0.8 0.8], [0 0 0], 'Time (s)', 'Speed (mm/s)', '', filename);
plot_tj_individuals(t, dzspeedwithdraw_I, [0 1 1], [0 0 1], 'Time (s)', 'Speed (mm/s)', '', filename);
yl = ylim;
plot([0 0], yl, '--k');
figure;
plot_tj_acrosstrial(t, dzspeedwithdraw_NI, [0 0 0], [0 0 0], 'Speed (mm/s)');
plot_tj_acrosstrial(t, dzspeedwithdraw_I, [0 0 1], [0 0 1], 'Speed (mm/s)');
yl = ylim;
plot([0 0], yl, '--k');
figure;
plot_tj_individuals(t, dzspeedwithdraw_diff, [0.8 0.8 0.8], [0 0 0], 'Time (s)', 'dSpeed (mm/s)', '', filename);
yl = ylim;
plot([0 0], yl, '--k');
plot([t(1) t(end)], [0 0], '--k');