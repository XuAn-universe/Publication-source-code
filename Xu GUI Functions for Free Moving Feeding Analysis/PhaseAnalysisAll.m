function PhaseAnalysisAll(app, Exp_Path, FrameRate, nmedian)
value = app.TrialsListBox.Value;
value = sort(value);

load([Exp_Path '\Analysis_Session.mat'], 'Video_annotation');

% process trajectories
if app.TrackingDataCheckBox.Value
    label_table = table2array(app.UITable.Data);
else
    helpdlg('''Tracking Data'' needs to be checked');
    return;
end
pointID = app.DropDown.Value;

try
    audiolocation = Exp_Path(1:end-7);
    temp = load([audiolocation '\Detected_Bite_Events.mat']);
    Bite_events = temp.Audio_analysis;
catch
    errordlg('Please detect the bites first.', 'Error');
    return;
end

if app.sOnTwiceButton.Value
    ninhibition = 2;
    tinhibition = [4 13; 8 17];
elseif app.NoLightButton.Value || app.WholeTrialOnButton.Value || app.BiteDeviceButton.Value
    ninhibition = 1;
    tinhibition = [0 inf];
elseif app.sDelay4sOnButton.Value
    ninhibition = 1;
    tinhibition = [4 8];
elseif app.sDelay4sOnButton_2.Value
    ninhibition = 1;
    tinhibition = [0 4];
end
FrameRate = round(FrameRate);
time_max = inf; % 12 or inf

phase_NI = [];
phase_I = [];

for i = value
    disp(['Processing trial ' num2str(i)]);
    events_related2drop = [];
    try
        temp = load([Exp_Path '\LabelledEvents' num2str(i) '.mat']);
        LabelledEvents = temp.LabelledEvents;
        MouthOpen = LabelledEvents.MouthOpen;
        TongueOut = LabelledEvents.TongueOut;
        PawLReachStart = LabelledEvents.PawLReachStart;
        PawRReachStart = LabelledEvents.PawRReachStart;
        events_related2drop = [MouthOpen; TongueOut; PawLReachStart; PawRReachStart];
    catch
        errordlg([Exp_Path '\LabelledEvents' num2str(i) '.mat is missing!'], 'Error');
    end
    
    if app.sDelay4sOnButton.Value
        if any(events_related2drop >= 4 & events_related2drop <= 8) % mouse drops pasta during inhibition
            continue;
        end
    end

    [y_PG1, y_PG3, z_PG1, z_PG2, z_PG3, x, y, z, speed, acceleration, laserstart, laserstop] =...
        trajectory_postprocessing(pointID, Exp_Path, i, label_table, nmedian, FrameRate);
    
    [y_anklel_PG1, y_anklel_PG3, z_anklel_PG1, z_anklel_PG2, z_anklel_PG3, x_anklel, y_anklel, z_anklel, speed_anklel, acceleration_anklel, laserstart, laserstop] =...
        trajectory_postprocessing(17, Exp_Path, i, label_table, nmedian, FrameRate);
    
    [y_ankler_PG1, y_ankler_PG3, z_ankler_PG1, z_ankler_PG2, z_ankler_PG3, x_ankler, y_ankler, z_ankler, speed_ankler, acceleration_ankler, laserstart, laserstop] =...
        trajectory_postprocessing(18, Exp_Path, i, label_table, nmedian, FrameRate);
    
    bite_timestamps = Bite_events(i).time_bites;
    if app.NoLightButton.Value
        laserstart = [];
        laserstop = [];
    else
        laser_timestamps = Bite_events(i).laser_timestamps;
    end
    
    time_ready2bite = Video_annotation(i).time_ready2bite;
    time_feeding_end = Video_annotation(i).time_feeding_end;
    if isempty(time_feeding_end)
        time_feeding_end = time_max;
    end
    time_range = [min(bite_timestamps) min(time_feeding_end, time_max)];
    
    t = (1:size(z, 1))'/FrameRate;
    time_rangeID = (t >= time_range(1) & t <= time_range(2));
    t_signal = t(time_rangeID);
    try
        if pointID == 10 % left paw
            signal = z(time_rangeID, 1)-z_anklel(time_rangeID, 1);
        elseif pointID == 16 % right paw
            signal = z(time_rangeID, 1)-z_ankler(time_rangeID, 1);
        end
        filter_signal = FIR_bandpass_filter(signal, FrameRate, [0.4 10], 0.1, nmedian, 0); % FIR_bandpass_filter(signal, SampleRate, passband, transition_width, nan_window, check_filter)
    catch
        if pointID == 10 % left paw
            signal = z(time_rangeID, 1)-z_ankler(time_rangeID, 1);
        elseif pointID == 16 % right paw
            signal = z(time_rangeID, 1)-z_anklel(time_rangeID, 1);
        end
        filter_signal = FIR_bandpass_filter(signal, FrameRate, [0.4 10], 0.1, nmedian, 0); % FIR_bandpass_filter(signal, SampleRate, passband, transition_width, nan_window, check_filter)
    end
    
    % compute phase
    hilbert_data = hilbert(filter_signal);
    phase_data = angle(hilbert_data)/pi*180; % this angle is the cosine angle
    phase_data(phase_data < 0) = phase_data(phase_data < 0)+360; % adjust phase to 0-360

    if ~isempty(laserstart)
        if ~isempty(bite_timestamps)
            for k = 1:ninhibition
                for m = 1:numel(bite_timestamps)
                    timestamp = bite_timestamps(m);
                    if timestamp >= laser_timestamps(2*k-1) && timestamp <= laser_timestamps(2*k)
                        [~, id] = min(abs(t_signal-timestamp));
                        phase_I = [phase_I phase_data(id)];
                    end
                end
            end
        end
    else
        if ~isempty(bite_timestamps)
            for k = 1:ninhibition
                for m = 1:numel(bite_timestamps)
                    timestamp = bite_timestamps(m);
                    if timestamp >= tinhibition(2*k-1) && timestamp <= tinhibition(2*k)
                        [~, id] = min(abs(t_signal-timestamp));
                        phase_NI = [phase_NI phase_data(id)];
                    end
                end
            end
        end
    end
end

result.phase_NI = phase_NI;
result.phase_I = phase_I;
assignin('base', 'result', result);

if ~app.NoLightButton.Value
    step = 9;
    edges = 0:step:360;
    N_NI = histcounts(phase_NI, edges, 'Normalization', 'probability');
    N_I = histcounts(phase_I, edges, 'Normalization', 'probability');
    figure;
    subplot(1, 2, 1);
    polarhistogram(circ_ang2rad(phase_NI), edges/180*pi, 'FaceColor', [0 0 0], 'EdgeColor', [0 0 0], 'FaceAlpha', 0.5, 'Normalization', 'probability');
    hold on;
    polarhistogram(circ_ang2rad(phase_I), edges/180*pi, 'FaceColor', [0 1 1], 'EdgeColor', [0 1 1], 'FaceAlpha', 0.5, 'Normalization', 'probability');
    set(gca, 'FontSize', 12, 'ButtonDownFcn', @extract_figure);
    
    subplot(1, 2, 2);
    histogram(phase_NI, edges, 'FaceColor', [0 0 0], 'EdgeColor', [0 0 0], 'FaceAlpha', 1, 'Normalization', 'probability');
    hold on;
    histogram(phase_I, edges, 'FaceColor', 'none', 'EdgeColor', [0 1 1], 'FaceAlpha', 1, 'Normalization', 'probability');
    plot(edges, max([N_NI N_I])/2*cosd(edges)+max([N_NI N_I])/2, '-r');
    xlim([0 360]);
    set(gca, 'TickDir', 'out', 'FontSize', 12, 'ButtonDownFcn', @extract_figure);
    box off;
    xlabel(['Phase (' char(176) ')']);
    ylabel('Probability');
    
    [pval, table] = circ_wwtest(circ_ang2rad(phase_NI), circ_ang2rad(phase_I));
    disp(['WW test: p = ' num2str(pval)]);
    [pval, med, P] = circ_cmtest(circ_ang2rad(phase_NI), circ_ang2rad(phase_I));
    disp(['CM test: p = ' num2str(pval)]);
    [pval, k, K] = circ_kuipertest(circ_ang2rad(phase_NI), circ_ang2rad(phase_I));
    disp(['Kuiper test: p < ' num2str(pval)]);
else
    step = 9;
    edges = 0:step:360;
    N_NI = histcounts(phase_NI, edges, 'Normalization', 'probability');
    figure;
    subplot(1, 2, 1);
    polarhistogram(circ_ang2rad(phase_NI), edges/180*pi, 'FaceColor', [0 0 0], 'EdgeColor', [0 0 0], 'FaceAlpha', 0.5, 'Normalization', 'probability');
    set(gca, 'FontSize', 12, 'ButtonDownFcn', @extract_figure);
    
    subplot(1, 2, 2);
    histogram(phase_NI, edges, 'FaceColor', [0 0 0], 'EdgeColor', [0 0 0], 'FaceAlpha', 1, 'Normalization', 'probability');
    hold on;
    plot(edges, max(N_NI)/2*cosd(edges)+max(N_NI)/2, '-r');
    xlim([0 360]);
    set(gca, 'TickDir', 'out', 'FontSize', 12, 'ButtonDownFcn', @extract_figure);
    box off;
    xlabel(['Phase (' char(176) ')']);
    ylabel('Probability');
end