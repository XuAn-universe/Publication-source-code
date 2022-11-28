function AHOrientation(app, Exp_Path, FrameRate, nmedian)
trial = app.TrialsListBox.Value;
if numel(trial) ~= 1
    helpdlg('You can only choose one trial');
    return;
end
% get bite events
audiolocation = Exp_Path(1:end-7);
try
    temp = load([audiolocation '\Detected_Bite_Events.mat']);
    Bite_events = temp.Audio_analysis;
    bite_timestamps = Bite_events(trial).time_bites;
    bite_amplitudes = Bite_events(trial).amplitude_bites/(max(Bite_events(trial).amplitude_bites));
catch
    bite_timestamps = [];
    bite_amplitudes = [];
end
% process trajectories
if app.TrackingDataCheckBox.Value
    label_table = table2array(app.UITable.Data);
else
    helpdlg('''Tracking Data'' needs to be checked');
    return;
end


% [~, ~, ~, ~, ~, x_top, y_top, z_top, ~, ~, laserstart, laserstop] = trajectory_postprocessing(21, Exp_Path, trial,...
%     label_table, nmedian, FrameRate);
% [~, ~, ~, ~, ~, x_bottom, y_bottom, z_bottom, ~, ~, ~, ~] = trajectory_postprocessing(22, Exp_Path, trial,...
%     label_table, nmedian, FrameRate);
% [~, ~, ~, ~, ~, x_center, y_center, z_center, ~, ~, ~, ~] = trajectory_postprocessing(33, Exp_Path, trial,...
%     label_table, nmedian, FrameRate);
% x = [x_top(:, 1) x_bottom(:, 1) x_center(:, 1)];
% y = [y_top(:, 1) y_bottom(:, 1) y_center(:, 1)];
% z = [z_top(:, 1) z_bottom(:, 1) z_center(:, 1)];
% orientationxy = nan(size(x));
% orientationxy(:, 1) = atand((z(:, 1)-z(:, 2))./(sqrt((x(:, 1)-x(:, 2)).^2+(y(:, 1)-y(:, 2)).^2)));
% orientationxy(:, 2) = atand((z(:, 1)-z(:, 3))./(sqrt((x(:, 1)-x(:, 3)).^2+(y(:, 1)-y(:, 3)).^2)));
% orientationxy(:, 3) = atand((z(:, 2)-z(:, 3))./(sqrt((x(:, 2)-x(:, 3)).^2+(y(:, 2)-y(:, 3)).^2)));
% orientationxy = abs(orientationxy);
% orientationxy = mean(orientationxy, 2, 'omitnan');

[y_topPG1, y_topPG3, z_topPG1, z_topPG2, z_topPG3, x_top, ~, ~, ~, ~, laserstart, laserstop] = trajectory_postprocessing(21, Exp_Path, trial,...
    label_table, nmedian, FrameRate);
[y_bottomPG1, y_bottomPG3, z_bottomPG1, z_bottomPG2, z_bottomPG3, x_bottom, ~, ~, ~, ~, ~, ~] = trajectory_postprocessing(22, Exp_Path, trial,...
    label_table, nmedian, FrameRate);
[y_centerPG1, y_centerPG3, z_centerPG1, z_centerPG2, z_centerPG3, x_center, ~, ~, ~, ~, ~, ~] = trajectory_postprocessing(33, Exp_Path, trial,...
    label_table, nmedian, FrameRate);
x_all = [x_top(:, 1) x_center(:, 1) x_bottom(:, 1)];
y_allPG1 = [y_topPG1(:, 1) y_centerPG1(:, 1) y_bottomPG1(:, 1)];
y_allPG3 = [y_topPG3(:, 1) y_centerPG3(:, 1) y_bottomPG3(:, 1)];
z_allPG1 = [z_topPG1(:, 1) z_centerPG1(:, 1) z_bottomPG1(:, 1)];
z_allPG2 = [z_topPG2(:, 1) z_centerPG2(:, 1) z_bottomPG2(:, 1)];
z_allPG3 = [z_topPG3(:, 1) z_centerPG3(:, 1) z_bottomPG3(:, 1)];
avc = slope_estimation(x_all, z_allPG2);
% avc = mean(avc, 2, 'omitnan');
bvcPG1 = slope_estimation(y_allPG1, z_allPG1);
bvcPG3 = slope_estimation(y_allPG3, z_allPG3);
bvc = mean([bvcPG1 bvcPG3], 2, 'omitnan');
orientationxy = acosd(1./sqrt(avc.^2+bvc.^2+1));
orientationxy = 90-orientationxy;

orientationxy_filtered = medfilt1(orientationxy, nmedian, 'omitnan', 'truncate');
t = (1:size(orientationxy, 1))'/FrameRate;
time_range = [app.timerangesEditField.Value app.toEditField.Value];
if app.NoLightButton.Value || ~app.OptSessionCheckBox.Value
    laserstart = [];
    laserstop = [];
else
    laserstart = laserstart/FrameRate;
    laserstop = laserstop/FrameRate;
end

% raw and filtered orientation with and without bite events
figure;
subplot(1, 2, 1);
plot_tj_singletrial(t, [orientationxy orientationxy_filtered], time_range, laserstart, laserstop,...
    [0 1 0], {'r', 'k'}, ['Orientation (' char(176) ')'], {'XY', 'XY filtered'});
subplot(1, 2, 2);
plot_tj_singletrial(t, orientationxy, time_range, laserstart, laserstop,...
    [0 1 0], {'k'}, ['Orientation (' char(176) ')'], 'XY');
if ~isempty(bite_timestamps)
    plot_biteevents(bite_timestamps, bite_amplitudes);
end