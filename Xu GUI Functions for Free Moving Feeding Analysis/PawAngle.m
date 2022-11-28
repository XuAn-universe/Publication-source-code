function PawAngle(app, Exp_Path, FrameRate, nmedian, pawsegment, xyzspeed_threshold)
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

[~, ~, ~, ~, ~, x_PL, y_PL, z_PL, ~, ~, laserstart, laserstop] = trajectory_postprocessing(10, Exp_Path, trial,...
    label_table, nmedian, FrameRate);
[~, ~, ~, ~, ~, x_PR, y_PR, z_PR, ~, ~, ~, ~] = trajectory_postprocessing(16, Exp_Path, trial,...
    label_table, nmedian, FrameRate);
vPL = [x_PL(:, 1) y_PL(:, 1) z_PL(:, 1)];
vPL = circshift(vPL, -pawsegment, 1)-vPL;
vPL(end-pawsegment+1:end, :) = nan;
vPR = [x_PR(:, 1) y_PR(:, 1) z_PR(:, 1)];
vPR = circshift(vPR, -pawsegment, 1)-vPR;
vPR(end-pawsegment+1:end, :) = nan;
vangle = atan2d(vecnorm(cross(vPL, vPR, 2), 2, 2), dot(vPL, vPR, 2));
vangle_filtered = medfilt1(vangle, nmedian, 'omitnan', 'truncate');

% ID_exclude = sqrt(sum(vPL.^2, 2)) < xyzspeed_threshold/FrameRate*pawsegment | sqrt(sum(vPR.^2, 2)) < xyzspeed_threshold/FrameRate*pawsegment;
% vangle(ID_exclude) = nan;
% vangle_filtered(ID_exclude) = nan;

t = (1:size(vangle, 1))'/FrameRate;
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
plot_tj_singletrial(t, [vangle vangle_filtered], time_range, laserstart, laserstop,...
    [0 1 0], {'r', 'k'}, ['Angle (' char(176) ')'], {'Angle', 'Angle filtered'});
subplot(1, 2, 2);
plot_tj_singletrial(t, vangle, time_range, laserstart, laserstop,...
    [0 1 0], {'k'}, ['Angle (' char(176) ')'], 'Angle');
if ~isempty(bite_timestamps)
    plot_biteevents(bite_timestamps, bite_amplitudes);
end