function DetectSit(app, Exp_Path, FrameRate, nmedian)
trial = app.TrialsListBox.Value;

load([Exp_Path '\Analysis_Session.mat'], 'Video_annotation');
time_ready2bite = Video_annotation(trial).time_ready2bite;
if isempty(time_ready2bite)
    helpdlg('Mark ready to bite first!');
    return;
end

% process trajectories
if app.TrackingDataCheckBox.Value
    label_table = table2array(app.UITable.Data);
else
    helpdlg('''Tracking Data'' needs to be checked');
    return;
end

[~, ~, ~, ~, ~, ~, ~, z_nose, ~, ~, ~, ~] = trajectory_postprocessing(3, Exp_Path, trial, label_table, nmedian, FrameRate);
[~, ~, ~, ~, ~, ~, y_pawl, z_pawl, ~, ~, ~, ~] = trajectory_postprocessing(10, Exp_Path, trial, label_table, nmedian, FrameRate);

t = (1:size(z_nose, 1))'/FrameRate;

win = 5;
z_nose_smooth = smoothdata(z_nose(:, 1), 'gaussian', win, 'omitnan');
index_sit_nose = get_index_sit(z_nose_smooth, t, time_ready2bite);
z_pawl_smooth = smoothdata(z_pawl(:, 1), 'gaussian', win, 'omitnan');
index_sit_pawl_z = get_index_sit(z_pawl_smooth, t, time_ready2bite);

y_pawl_smooth = smoothdata(y_pawl(:, 1), 'gaussian', win, 'omitnan');
yspeed = [NaN; diff(y_pawl_smooth)];
yspeed_ready2bite = yspeed(t <= time_ready2bite);
y_pawl_smooth_ready2bite = y_pawl_smooth(t <= time_ready2bite);
[~, index_ymin] = min(y_pawl_smooth_ready2bite);
[~, fallingindex, ~, ~] = ZeroCrossingDetection(yspeed_ready2bite);
index_sit_pawl_y = fallingindex(find(fallingindex <= index_ymin, 1, 'last'));
index_sit_pawl = min(index_sit_pawl_y, index_sit_pawl_z);

figure;
subplot(1, 3, 1);
hold on;
plot(t, z_nose_smooth, '-k');
plot(t(index_sit_nose), z_nose_smooth(index_sit_nose), 'or');
xlabel('Time (s)');
ylabel('Z (mm)');
title('Nose');
set(gca, 'TickLength', [0 0], 'FontSize', 12, 'ButtonDownFcn', @extract_figure);

subplot(1, 3, 2);
hold on;
plot(t, z_pawl_smooth, '-k');
plot(t(index_sit_pawl_z), z_pawl_smooth(index_sit_pawl_z), 'sg');
plot(t(index_sit_pawl), z_pawl_smooth(index_sit_pawl), 'or');
xlabel('Time (s)');
ylabel('Z (mm)');
title('PawL');
set(gca, 'TickLength', [0 0], 'FontSize', 12, 'ButtonDownFcn', @extract_figure);

subplot(1, 3, 3);
hold on;
plot(t, y_pawl_smooth, '-k');
plot(t(index_sit_pawl_y), y_pawl_smooth(index_sit_pawl_y), 'sg');
plot(t(index_sit_pawl), y_pawl_smooth(index_sit_pawl), 'or');
xlabel('Time (s)');
ylabel('Y (mm)');
title('PawL');
set(gca, 'TickLength', [0 0], 'FontSize', 12, 'ButtonDownFcn', @extract_figure);

app.NoseEditField.Value = t(index_sit_nose);
app.PawLEditField.Value = t(index_sit_pawl);

function index_sit = get_index_sit(z, t, time_ready2bite)
tskip = 0.5;
zspeed = [NaN; diff(z)];
zspeed = zspeed(t <= time_ready2bite);
[risingindex, ~, ~, ~] = ZeroCrossingDetection(zspeed);
zrange = zeros(1, numel(risingindex));
for i = 1:numel(risingindex)
    if i ~= numel(risingindex)
        zrange(i) = range(z(risingindex(i):risingindex(i+1)));
    else
        zrange(i) = range(z(risingindex(i):find(t <= time_ready2bite, 1, 'last')));
    end
end
zrange(t(risingindex) <= tskip) = 0;
index_sit = risingindex(zrange == max(zrange));