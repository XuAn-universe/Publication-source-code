function PhaseAnalysis(app, Exp_Path, FrameRate, nmedian)
trial = app.TrialsListBox.Value;
if numel(trial) ~= 1
    helpdlg('You can only choose one trial');
    return;
end

FrameRate = round(FrameRate);

load([Exp_Path '\Analysis_Session.mat'], 'Video_annotation');
time_ready2bite = Video_annotation(trial).time_ready2bite;
time_feeding_end = Video_annotation(trial).time_feeding_end;
time_range = [time_ready2bite time_feeding_end];

% get bite events
audiolocation = Exp_Path(1:end-7);
try
    temp = load([audiolocation '\Detected_Bite_Events.mat']);
    Bite_events = temp.Audio_analysis;
    bite_timestamps = Bite_events(trial).time_bites;
%     bite_amplitudes = Bite_events(trial).amplitude_bites/(max(Bite_events(trial).amplitude_bites));
    bite_amplitudes = Bite_events(trial).amplitude_bites./(Bite_events(trial).amplitude_bites);
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
pointID = app.DropDown.Value;

[y_PG1, y_PG3, z_PG1, z_PG2, z_PG3, x, y, z, speed, acceleration, laserstart, laserstop] =...
    trajectory_postprocessing(pointID, Exp_Path, trial, label_table, nmedian, FrameRate);

[y_anklel_PG1, y_anklel_PG3, z_anklel_PG1, z_anklel_PG2, z_anklel_PG3, x_anklel, y_anklel, z_anklel, speed_anklel, acceleration_anklel, laserstart, laserstop] =...
    trajectory_postprocessing(17, Exp_Path, trial, label_table, nmedian, FrameRate);

[y_ankler_PG1, y_ankler_PG3, z_ankler_PG1, z_ankler_PG2, z_ankler_PG3, x_ankler, y_ankler, z_ankler, speed_ankler, acceleration_ankler, laserstart, laserstop] =...
    trajectory_postprocessing(18, Exp_Path, trial, label_table, nmedian, FrameRate);

t = (1:numel(y_PG1))'/FrameRate;
time_rangeID = (t >= time_range(1) & t <= time_range(2));

if app.NoLightButton.Value || ~app.OptSessionCheckBox.Value
    laserstart = [];
    laserstop = [];
else
    laserstart = laserstart/FrameRate;
    laserstop = laserstop/FrameRate;
end

% raw X, Y, Z trajectories with bite events
figure('Name', 'REF: left ankle');
subplot(1, 4, 1);
plot_tj_singletrial(t, z(:, 1)-z_anklel(:, 1), time_range, laserstart, laserstop, [0 1 0], {'k'}, 'dZ (mm)', {'dZ'});
if ~isempty(bite_timestamps)
    plot_biteevents(bite_timestamps, bite_amplitudes);
end

subplot(1, 4, 2);
plot_tj_singletrial(t, vecnorm([x(:, 1)-x_anklel(:, 1) z(:, 1)-z_anklel(:, 1)], 2, 2), time_range, laserstart, laserstop, [0 1 0], {'k'}, 'dXZ (mm)', {'dXZ'});
if ~isempty(bite_timestamps)
    plot_biteevents(bite_timestamps, bite_amplitudes);
end

subplot(1, 4, 3);
plot_tj_singletrial(t, vecnorm([y(:, 1)-y_anklel(:, 1) z(:, 1)-z_anklel(:, 1)], 2, 2), time_range, laserstart, laserstop, [0 1 0], {'k'}, 'dYZ (mm)', {'dYZ'});
if ~isempty(bite_timestamps)
    plot_biteevents(bite_timestamps, bite_amplitudes);
end

subplot(1, 4, 4);
plot_tj_singletrial(t, vecnorm([x(:, 1)-x_anklel(:, 1) y(:, 1)-y_anklel(:, 1) z(:, 1)-z_anklel(:, 1)], 2, 2), time_range, laserstart, laserstop, [0 1 0], {'k'}, 'dXYZ (mm)', {'dXYZ'});
if ~isempty(bite_timestamps)
    plot_biteevents(bite_timestamps, bite_amplitudes);
end

figure('Name', 'REF: right ankle');
subplot(1, 4, 1);
plot_tj_singletrial(t, z(:, 1)-z_ankler(:, 1), time_range, laserstart, laserstop, [0 1 0], {'k'}, 'dZ (mm)', {'dZ'});
if ~isempty(bite_timestamps)
    plot_biteevents(bite_timestamps, bite_amplitudes);
end

subplot(1, 4, 2);
plot_tj_singletrial(t, vecnorm([x(:, 1)-x_ankler(:, 1) z(:, 1)-z_ankler(:, 1)], 2, 2), time_range, laserstart, laserstop, [0 1 0], {'k'}, 'dXZ (mm)', {'dXZ'});
if ~isempty(bite_timestamps)
    plot_biteevents(bite_timestamps, bite_amplitudes);
end

subplot(1, 4, 3);
plot_tj_singletrial(t, vecnorm([y(:, 1)-y_ankler(:, 1) z(:, 1)-z_ankler(:, 1)], 2, 2), time_range, laserstart, laserstop, [0 1 0], {'k'}, 'dYZ (mm)', {'dYZ'});
if ~isempty(bite_timestamps)
    plot_biteevents(bite_timestamps, bite_amplitudes);
end

subplot(1, 4, 4);
plot_tj_singletrial(t, vecnorm([x(:, 1)-x_ankler(:, 1) y(:, 1)-y_ankler(:, 1) z(:, 1)-z_ankler(:, 1)], 2, 2), time_range, laserstart, laserstop, [0 1 0], {'k'}, 'dXYZ (mm)', {'dXYZ'});
if ~isempty(bite_timestamps)
    plot_biteevents(bite_timestamps, bite_amplitudes);
end

t_signal = t(time_rangeID);
ylabel_text = 'Signal (mm)';
signal_all = nan(sum(time_rangeID), 8);
signal_all(:, 1) = z(time_rangeID, 1)-z_anklel(time_rangeID, 1);
signal_all(:, 2) = vecnorm([x(time_rangeID, 1)-x_anklel(time_rangeID, 1) z(time_rangeID, 1)-z_anklel(time_rangeID, 1)], 2, 2);
signal_all(:, 3) = vecnorm([y(time_rangeID, 1)-y_anklel(time_rangeID, 1) z(time_rangeID, 1)-z_anklel(time_rangeID, 1)], 2, 2);
signal_all(:, 4) = vecnorm([x(time_rangeID, 1)-x_anklel(time_rangeID, 1) y(time_rangeID, 1)-y_anklel(time_rangeID, 1) z(time_rangeID, 1)-z_anklel(time_rangeID, 1)], 2, 2);
signal_all(:, 5) = z(time_rangeID, 1)-z_ankler(time_rangeID, 1);
signal_all(:, 6) = vecnorm([x(time_rangeID, 1)-x_ankler(time_rangeID, 1) z(time_rangeID, 1)-z_ankler(time_rangeID, 1)], 2, 2);
signal_all(:, 7) = vecnorm([y(time_rangeID, 1)-y_ankler(time_rangeID, 1) z(time_rangeID, 1)-z_ankler(time_rangeID, 1)], 2, 2);
signal_all(:, 8) = vecnorm([x(time_rangeID, 1)-x_ankler(time_rangeID, 1) y(time_rangeID, 1)-y_ankler(time_rangeID, 1) z(time_rangeID, 1)-z_ankler(time_rangeID, 1)], 2, 2);
fname = {'Z left ankle', 'XZ left ankle', 'YZ left ankle', 'XYZ left ankle', 'Z right ankle', 'XZ right ankle', 'YZ right ankle', 'XYZ right ankle'};
for i = 1:8
    try
        signal = signal_all(:, i);
        figure('Name', fname{i});
        subplot(1, 2, 1);
        filter_signal = FIR_bandpass_filter(signal, FrameRate, [0.4 10], 0.1, nmedian, 1); % FIR_bandpass_filter(signal, SampleRate, passband, transition_width, nan_window, check_filter)
        
        subplot(1, 2, 2);
        hp(1) = plot(t_signal, filter_signal, '-b');
        hold on;
        if ~isempty(bite_amplitudes)
            plot_biteevents(bite_timestamps, bite_amplitudes);
        end
        hp(2) = plot(t_signal, signal, '-k');
        plot([time_range(1) time_range(end)], [0 0], '-k');
        xlim(time_range);
        xlabel('Time (s)');
        ylabel(ylabel_text);
        title('From ready-to-bite to feeding end');
        set(gca, 'FontSize', 12, 'ButtonDownFcn', @extract_figure);
        box off;
        yl = ylim;
        if ~isempty(laserstart)
            for j = 1:numel(laserstart)
                line([laserstart(j) laserstart(j)], yl, 'Color', [0 1 0], 'LineStyle', '--', 'LineWidth', 1);
                line([laserstop(j) laserstop(j)], yl, 'Color', [0 1 0], 'LineStyle', '--', 'LineWidth', 1);
            end
        end
        legend(hp, {'Band-pass filtered', 'Original'}, 'Location', 'northeast', 'FontSize', 12);
        
        % compute phase
        hilbert_data = hilbert(filter_signal);
        phase_data = angle(hilbert_data)/pi*180; % this angle is the cosine angle
        phase_data(phase_data < 0) = phase_data(phase_data < 0)+360; % adjust phase to 0-360
        
        figure('Name', fname{i});
        subplot(2, 2, 1);
        plot(t_signal, phase_data, '-k');
        hold on;
        phase_bite = plot_biteevents(bite_timestamps, bite_amplitudes); % this also computes bite phase
        xlim(time_range);
        xlabel('Time (s)');
        ylabel(['Phase (' char(176) ')']);
        title('From ready-to-bite to feeding end');
        set(gca, 'FontSize', 12, 'ButtonDownFcn', @extract_figure);
        box off;
        yl = ylim;
        if ~isempty(laserstart)
            for j = 1:numel(laserstart)
                line([laserstart(j) laserstart(j)], yl, 'Color', [0 1 0], 'LineStyle', '--', 'LineWidth', 1);
                line([laserstop(j) laserstop(j)], yl, 'Color', [0 1 0], 'LineStyle', '--', 'LineWidth', 1);
            end
        end
        
        subplot(2, 2, 2);
        scatter(bite_amplitudes, phase_bite, 36, [0 0 0], 'o');
        xlabel('Bite amplitude');
        ylabel(['Phase (' char(176) ')']);
        set(gca, 'FontSize', 12, 'ButtonDownFcn', @extract_figure);
        box off;
        
        step = 9;
        edges = 0:step:360;
        N = histcounts(phase_bite, edges, 'Normalization', 'probability');
        subplot(2, 2, 3);
        polarhistogram(circ_ang2rad(phase_bite), edges/180*pi, 'FaceColor', [0 0 0], 'EdgeColor', [0 0 0], 'FaceAlpha', 0.5, 'Normalization', 'probability');
        set(gca, 'FontSize', 12, 'ButtonDownFcn', @extract_figure);
        
        subplot(2, 2, 4);
        histogram(phase_bite, edges, 'FaceColor', [0 0 0], 'FaceAlpha', 1, 'Normalization', 'probability');
        hold on;
        plot(edges, max(N)/2*cosd(edges)+max(N)/2, '-c');
        xlim([0 360]);
        set(gca, 'TickDir', 'out', 'FontSize', 12, 'ButtonDownFcn', @extract_figure);
        box off;
        xlabel(['Phase (' char(176) ')']);
        ylabel('Probability');
        
        [mu, ul, ll] = circ_mean(circ_ang2rad(phase_bite'));
        r = circ_r(circ_ang2rad(phase_bite'));
        title(['Mean = ' num2str(mu/pi*180) '; r = ' num2str(r)]);
    end
end