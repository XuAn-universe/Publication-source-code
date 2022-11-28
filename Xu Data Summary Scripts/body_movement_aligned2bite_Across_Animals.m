%*---------------------------------------------------------------------*
% Huang Lab
% Duke University
% Author : Xu An, May 2022
% xu.an@duke.edu
% Version: 1.0
%*---------------------------------------------------------------------*
%% body movement aligned to bite
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
trange = 0.5;
trange = trange*FrameRate;
t = (-trange:1:trange)'/FrameRate;

N = numel(filename);
zbite_NI = nan(2*trange+1, N);
zbite_I = nan(2*trange+1, N);
zspeedbite_NI = nan(2*trange+1, N);
zspeedbite_I = nan(2*trange+1, N);
for i = 1:N
    temp = load([pathname filename{i}]);
    result = temp.result;
    zbite_NI(:, i) = mean(result.zbite_NI{1}, 2, 'omitnan');
    if isfield(result, 'zbite_I') && ~isempty(result.zbite_I{1})
        zbite_I(:, i) = mean(result.zbite_I{1}, 2, 'omitnan');
    end
    zspeedbite_NI(:, i) = mean(result.zspeedbite_NI{1}, 2, 'omitnan');
    if isfield(result, 'zspeedbite_I') && ~isempty(result.zspeedbite_I{1})
        zspeedbite_I(:, i) = mean(result.zspeedbite_I{1}, 2, 'omitnan');
    end
end
zbite_diff = zbite_I-zbite_NI;
% zbite_NI = zbite_NI-ones(2*trange+1, 1)*zbite_NI(trange+1, :); % normalization
% zbite_I = zbite_I-ones(2*trange+1, 1)*zbite_I(trange+1, :); % normalization
figure;
plot_tj_individuals(t, zbite_NI, [0.8 0.8 0.8], [0 0 0], 'Time (s)', 'Z (mm)', '', filename);
plot_tj_individuals(t, zbite_I, [0 1 1], [0 0 1], 'Time (s)', 'Z (mm)', '', filename);
yl = ylim;
plot([0 0], yl, '--k');
figure;
plot_tj_acrosstrial(t, zbite_NI, [0 0 0], [0 0 0], 'Z (mm)');
plot_tj_acrosstrial(t, zbite_I, [0 0 1], [0 0 1], 'Z (mm)');
yl = ylim;
plot([0 0], yl, '--k');
figure;
plot_tj_individuals(t, zbite_diff, [0.8 0.8 0.8], [0 0 0], 'Time (s)', 'dZ (mm)', '', filename);
yl = ylim;
plot([0 0], yl, '--k');
plot([t(1) t(end)], [0 0], '--k');

zspeedbite_diff = zspeedbite_I-zspeedbite_NI;
figure;
plot_tj_individuals(t, zspeedbite_NI, [0.8 0.8 0.8], [0 0 0], 'Time (s)', 'Speed (mm/s)', '', filename);
plot_tj_individuals(t, zspeedbite_I, [0 1 1], [0 0 1], 'Time (s)', 'Speed (mm/s)', '', filename);
yl = ylim;
plot([0 0], yl, '--k');
figure;
plot_tj_acrosstrial(t, zspeedbite_NI, [0 0 0], [0 0 0], 'Speed (mm/s)');
plot_tj_acrosstrial(t, zspeedbite_I, [0 0 1], [0 0 1], 'Speed (mm/s)');
yl = ylim;
plot([0 0], yl, '--k');
figure;
plot_tj_individuals(t, zspeedbite_diff, [0.8 0.8 0.8], [0 0 0], 'Time (s)', 'dSpeed (mm/s)', '', filename);
yl = ylim;
plot([0 0], yl, '--k');
plot([t(1) t(end)], [0 0], '--k');

%% combine body parts together
figure;
plot_tj_acrosstrial(t, zbite_NI_nose-ones(2*trange+1, 1)*zbite_NI_nose(trange+1, :), [1 0 0], [1 0 0], 'Z (mm)');
plot_tj_acrosstrial(t, zbite_NI_gh-ones(2*trange+1, 1)*zbite_NI_nose(trange+1, :), [0 0 1], [0 0 1], 'Z (mm)');
plot_tj_acrosstrial(t, zbite_NI_sh-ones(2*trange+1, 1)*zbite_NI_nose(trange+1, :), [0 1 0], [0 1 0], 'Z (mm)');
yl = ylim;
plot([0 0], yl, '--k');

figure;
plot_tj_acrosstrial(t, zspeedbite_NI_nose, [1 0 0], [1 0 0], 'Speed (mm/s)');
plot_tj_acrosstrial(t, zspeedbite_NI_gh, [0 0 1], [0 0 1], 'Speed (mm/s)');
plot_tj_acrosstrial(t, zspeedbite_NI_sh, [0 1 0], [0 1 0], 'Speed (mm/s)');
yl = ylim;
plot([0 0], yl, '--k');

%% relative body movement aligned to bite
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
trange = 0.5;
trange = trange*FrameRate;
t = (-trange:1:trange)'/FrameRate;

N = numel(filename);
dzbite_NI = nan(2*trange+1, N);
dzbite_I = nan(2*trange+1, N);
orientationxybite_NI = nan(2*trange+1, N);
orientationxybite_I = nan(2*trange+1, N);
for i = 1:N
    temp = load([pathname filename{i}]);
    result = temp.result;
    dzbite_NI(:, i) = mean(result.dzbite_NI{1}, 2, 'omitnan');
    if isfield(result, 'dzbite_I') && ~isempty(result.dzbite_I{1})
        dzbite_I(:, i) = mean(result.dzbite_I{1}, 2, 'omitnan');
    end
    orientationxybite_NI(:, i) = mean(result.orientationxybite_NI{1}, 2, 'omitnan');
    if isfield(result, 'orientationxybite_I') && ~isempty(result.orientationxybite_I{1})
        orientationxybite_I(:, i) = mean(result.orientationxybite_I{1}, 2, 'omitnan');
    end
end
dzbite_diff = dzbite_I-dzbite_NI;
% dzbite_NI = dzbite_NI-ones(2*trange+1, 1)*dzbite_NI(trange+1, :); % normalization
% dzbite_I = dzbite_I-ones(2*trange+1, 1)*dzbite_I(trange+1, :); % normalization
figure;
plot_tj_individuals(t, dzbite_NI, [0.8 0.8 0.8], [0 0 0], 'Time (s)', 'dZ (mm)', '', filename);
plot_tj_individuals(t, dzbite_I, [0 1 1], [0 0 1], 'Time (s)', 'dZ (mm)', '', filename);
yl = ylim;
plot([0 0], yl, '--k');
figure;
plot_tj_acrosstrial(t, dzbite_NI, [0 0 0], [0 0 0], 'dZ (mm)');
plot_tj_acrosstrial(t, dzbite_I, [0 0 1], [0 0 1], 'dZ (mm)');
yl = ylim;
plot([0 0], yl, '--k');
figure;
plot_tj_individuals(t, dzbite_diff, [0.8 0.8 0.8], [0 0 0], 'Time (s)', 'dZ (mm)', '', filename);
yl = ylim;
plot([0 0], yl, '--k');
plot([t(1) t(end)], [0 0], '--k');

orientationxybite_diff = orientationxybite_I-orientationxybite_NI;
% orientationxybite_NI = orientationxybite_NI-ones(2*trange+1, 1)*orientationxybite_NI(trange+1, :); % normalization
% orientationxybite_I = orientationxybite_I-ones(2*trange+1, 1)*orientationxybite_I(trange+1, :); % normalization
figure;
plot_tj_individuals(t, orientationxybite_NI, [0.8 0.8 0.8], [0 0 0], 'Time (s)', ['Orientation (' char(176) ')'], '', filename);
plot_tj_individuals(t, orientationxybite_I, [0 1 1], [0 0 1], 'Time (s)', ['Orientation (' char(176) ')'], '', filename);
yl = ylim;
plot([0 0], yl, '--k');
figure;
plot_tj_acrosstrial(t, orientationxybite_NI, [0 0 0], [0 0 0], ['Orientation (' char(176) ')']);
plot_tj_acrosstrial(t, orientationxybite_I, [0 0 1], [0 0 1], ['Orientation (' char(176) ')']);
yl = ylim;
plot([0 0], yl, '--k');
figure;
plot_tj_individuals(t, orientationxybite_diff, [0.8 0.8 0.8], [0 0 0], 'Time (s)', ['dOrientation (' char(176) ')'], '', filename);
yl = ylim;
plot([0 0], yl, '--k');
plot([t(1) t(end)], [0 0], '--k');

%% combine two hands together
figure;
plot_tj_acrosstrial(t, dzbite_NI_gh-ones(2*trange+1, 1)*dzbite_NI_gh(trange+1, :), [0 0 1], [0 0 1], 'dZ (mm)');
plot_tj_acrosstrial(t, dzbite_NI_sh-ones(2*trange+1, 1)*dzbite_NI_sh(trange+1, :), [0 1 0], [0 1 0], 'dZ (mm)');
yl = ylim;
plot([0 0], yl, '--k');

figure;
plot_tj_acrosstrial(t, orientationxybite_NI_gh-ones(2*trange+1, 1)*orientationxybite_NI_gh(trange+1, :), [0 0 1], [0 0 1], ['Orientation (' char(176) ')']);
plot_tj_acrosstrial(t, orientationxybite_NI_sh-ones(2*trange+1, 1)*orientationxybite_NI_sh(trange+1, :), [0 1 0], [0 1 0], ['Orientation (' char(176) ')']);
yl = ylim;
plot([0 0], yl, '--k');