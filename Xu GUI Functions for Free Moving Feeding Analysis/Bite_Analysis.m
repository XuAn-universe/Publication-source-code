function Bite_Analysis(app, Exp_Path, bite_bin)
audiolocation = Exp_Path(1:end-7);
try
    temp = load([audiolocation '\Detected_Bite_Events.mat']);
    Bite_events = temp.Audio_analysis;
    temp = load([Exp_Path '\Analysis_Session.mat']);
    video_annotation = temp.Video_annotation;
    
    value = app.ListBox.Value;
    if ~isempty(value)
        for i = 1:numel(value)
            audiolocation = app.ListBox.Items{value(i)}(1:end-7);
            temp = load([audiolocation '\Detected_Bite_Events.mat']);
            Bite_events = [Bite_events; temp.Audio_analysis];
            
            temp = load([app.ListBox.Items{value(i)} '\Analysis_Session.mat']);
            if ~isfield(temp.Video_annotation, 'Omission')
                [temp.Video_annotation(:).Omission] = deal(0);
            end
            video_annotation = [video_annotation; temp.Video_annotation];
        end
    end
catch
    errordlg('Please detect the bites first.', 'Error');
    return;
end
if app.sDelay4sOnButton.Value
    firstafter2lastbefore_I = [];
    firstafter2lastbefore_NI = [];
elseif app.sDelay4sOnButton_2.Value || app.WholeTrialOnButton.Value || app.BiteDeviceButton.Value
    firstbite_time_I = [];
    firstbite_time_NI = [];
    if app.BiteDeviceButton.Value
        bite_duration_I = [];
        bite_duration_NI = [];
    end
elseif app.NoLightButton.Value
    firstbite_time = [];
    all_trials = [];
    totalbites = [];
    bite_intervals = [];
    bite_amplitudes = [];
else
    firstafter2lastbefore_I1 = [];
    firstafter2lastbefore_NI1 = [];
    firstafter2lastbefore_I2 = [];
    firstafter2lastbefore_NI2 = [];
    totalbites_I1 = [];
    totalbites_NI1 = [];
    bite_intervals_I1 = [];
    bite_intervals_NI1 = [];
    bite_amplitudes_I1 = [];
    bite_amplitudes_NI1 = [];
    totalbites_I2 = [];
    totalbites_NI2 = [];
    bite_intervals_I2 = [];
    bite_intervals_NI2 = [];
    bite_amplitudes_I2 = [];
    bite_amplitudes_NI2 = [];
    trials_I = [];
    trials_NI = [];
end
if app.sDelay4sOnButton.Value || app.sDelay4sOnButton_2.Value || app.WholeTrialOnButton.Value || app.BiteDeviceButton.Value
    totalbites_I = [];
    totalbites_NI = [];
    bite_intervals_I = [];
    bite_intervals_NI = [];
    bite_amplitudes_I = [];
    bite_amplitudes_NI = [];
    trials_I = [];
    trials_NI = [];
end

max_timestamp = 0;
for i = 1:numel(video_annotation)
    if ~video_annotation(i).Disgard
        bite_timestamps = Bite_events(i).time_bites;
        bite_amplitudes = Bite_events(i).amplitude_bites/(max(Bite_events(i).amplitude_bites));
        if ~isempty(bite_timestamps)
            max_timestamp = max(max_timestamp, max(bite_timestamps));
        end
        if app.sDelay4sOnButton.Value
            if isempty(Bite_events(i).laser_timestamps)
                if isempty(trials_NI)
                    trials_NI(1).bite_timestamps = bite_timestamps;
                    trials_NI(1).bite_amplitudes = bite_amplitudes;
                    if ~isempty(bite_timestamps)
                        firstafter2lastbefore_NI(1) = bite_timestamps(find(bite_timestamps >= 4, 1))-...
                            bite_timestamps(find(bite_timestamps < 4, 1, 'last'));
                    end
                    totalbites_NI(1) = numel(find(bite_timestamps >= 4 & bite_timestamps <= 8));
                else
                    trials_NI(end+1).bite_timestamps = bite_timestamps;
                    trials_NI(end).bite_amplitudes = bite_amplitudes;
                    if ~isempty(bite_timestamps)
                        firstafter2lastbefore_NI(end+1) = bite_timestamps(find(bite_timestamps >= 4, 1))-...
                            bite_timestamps(find(bite_timestamps < 4, 1, 'last'));
                    end
                    totalbites_NI(end+1) = numel(find(bite_timestamps >= 4 & bite_timestamps <= 8));
                end
                bite_intervals_NI = [bite_intervals_NI diff(bite_timestamps(bite_timestamps >= 4 & bite_timestamps <= 8))];
                bite_amplitudes_NI = [bite_amplitudes_NI bite_amplitudes(bite_timestamps >= 4 & bite_timestamps <= 8)];
            else
                if isempty(trials_I)
                    trials_I(1).bite_timestamps = bite_timestamps;
                    trials_I(1).bite_amplitudes = bite_amplitudes;
                    trials_I(1).laser_timestamps = Bite_events(i).laser_timestamps;
                    if ~isempty(bite_timestamps)
                        firstafter2lastbefore_I(1) = bite_timestamps(find(bite_timestamps >= Bite_events(i).laser_timestamps(1), 1))-...
                            bite_timestamps(find(bite_timestamps < Bite_events(i).laser_timestamps(1), 1, 'last'));
                    end
                    totalbites_I(1) = numel(find(bite_timestamps >= Bite_events(i).laser_timestamps(1) &...
                        bite_timestamps <= Bite_events(i).laser_timestamps(2)));
                else
                    trials_I(end+1).bite_timestamps = bite_timestamps;
                    trials_I(end).bite_amplitudes = bite_amplitudes;
                    trials_I(end).laser_timestamps = Bite_events(i).laser_timestamps;
                    if ~isempty(bite_timestamps)
                        firstafter2lastbefore_I(end+1) = bite_timestamps(find(bite_timestamps >= Bite_events(i).laser_timestamps(1), 1))-...
                            bite_timestamps(find(bite_timestamps < Bite_events(i).laser_timestamps(1), 1, 'last'));
                    end
                    totalbites_I(end+1) = numel(find(bite_timestamps >= Bite_events(i).laser_timestamps(1) &...
                        bite_timestamps <= Bite_events(i).laser_timestamps(2)));
                end
                bite_intervals_I = [bite_intervals_I diff(bite_timestamps(bite_timestamps >= Bite_events(i).laser_timestamps(1) &...
                    bite_timestamps <= Bite_events(i).laser_timestamps(2)))];
                bite_amplitudes_I = [bite_amplitudes_I bite_amplitudes(bite_timestamps >= Bite_events(i).laser_timestamps(1) &...
                    bite_timestamps <= Bite_events(i).laser_timestamps(2))];
            end
        elseif app.sDelay4sOnButton_2.Value
            if isempty(Bite_events(i).laser_timestamps)
                if isempty(trials_NI)
                    trials_NI(1).bite_timestamps = bite_timestamps;
                    trials_NI(1).bite_amplitudes = bite_amplitudes;
                    if ~isempty(bite_timestamps)
                        firstbite_time_NI(1) = bite_timestamps(1);
                    end
                    totalbites_NI(1) = numel(find(bite_timestamps >= 0 & bite_timestamps <= 4));
                else
                    trials_NI(end+1).bite_timestamps = bite_timestamps;
                    trials_NI(end).bite_amplitudes = bite_amplitudes;
                    if ~isempty(bite_timestamps)
                        firstbite_time_NI(end+1) = bite_timestamps(1);
                    end
                    totalbites_NI(end+1) = numel(find(bite_timestamps >= 0 & bite_timestamps <= 4));
                end
                bite_intervals_NI = [bite_intervals_NI diff(bite_timestamps(bite_timestamps >= 0 & bite_timestamps <= 4))];
                bite_amplitudes_NI = [bite_amplitudes_NI bite_amplitudes(bite_timestamps >= 0 & bite_timestamps <= 4)];
            else
                if isempty(trials_I)
                    trials_I(1).bite_timestamps = bite_timestamps;
                    trials_I(1).bite_amplitudes = bite_amplitudes;
                    trials_I(1).laser_timestamps = Bite_events(i).laser_timestamps;
                    if ~isempty(bite_timestamps)
                        firstbite_time_I(1) = bite_timestamps(1);
                    end
                    totalbites_I(1) = numel(find(bite_timestamps >= Bite_events(i).laser_timestamps(1) &...
                        bite_timestamps <= Bite_events(i).laser_timestamps(2)));
                else
                    trials_I(end+1).bite_timestamps = bite_timestamps;
                    trials_I(end).bite_amplitudes = bite_amplitudes;
                    trials_I(end).laser_timestamps = Bite_events(i).laser_timestamps;
                    if ~isempty(bite_timestamps)
                        firstbite_time_I(end+1) = bite_timestamps(1);
                    end
                    totalbites_I(end+1) = numel(find(bite_timestamps >= Bite_events(i).laser_timestamps(1) &...
                        bite_timestamps <= Bite_events(i).laser_timestamps(2)));
                end
                bite_intervals_I = [bite_intervals_I diff(bite_timestamps(bite_timestamps >= Bite_events(i).laser_timestamps(1) &...
                    bite_timestamps <= Bite_events(i).laser_timestamps(2)))];
                bite_amplitudes_I = [bite_amplitudes_I bite_amplitudes(bite_timestamps >= Bite_events(i).laser_timestamps(1) &...
                    bite_timestamps <= Bite_events(i).laser_timestamps(2))];
            end
        elseif app.NoLightButton.Value
            if isempty(all_trials)
                if ~isempty(bite_timestamps)
                    firstbite_time(1) = bite_timestamps(1);
                end
                all_trials(1).bite_timestamps = bite_timestamps;
                all_trials(1).bite_amplitudes = bite_amplitudes;
                totalbites(1) = numel(bite_timestamps);
            else
                if ~isempty(bite_timestamps)
                    firstbite_time(end+1) = bite_timestamps(1);
                end
                all_trials(end+1).bite_timestamps = bite_timestamps;
                all_trials(end).bite_amplitudes = bite_amplitudes;
                totalbites(end+1) = numel(bite_timestamps);
            end
            bite_intervals = [bite_intervals diff(bite_timestamps)];
            bite_amplitudes = [bite_amplitudes bite_amplitudes];
        elseif app.WholeTrialOnButton.Value || app.BiteDeviceButton.Value
            if isempty(Bite_events(i).laser_timestamps)
                if isempty(trials_NI)
                    trials_NI(1).bite_timestamps = bite_timestamps;
                    trials_NI(1).bite_amplitudes = bite_amplitudes;
                    if ~isempty(bite_timestamps)
                        firstbite_time_NI(1) = bite_timestamps(1);
                    end
                    totalbites_NI(1) = numel(bite_timestamps);
                    if app.BiteDeviceButton.Value
                        bite_duration_NI(1) = 0;
                        for j = 1:size(Bite_events(i).bite_window, 1)
                            bite_timestamps_temp = bite_timestamps(bite_timestamps > Bite_events(i).bite_window(j, 1) & bite_timestamps < Bite_events(i).bite_window(j, 2));
                            bite_duration_NI(1) = bite_duration_NI(1)+range(bite_timestamps_temp);
                            bite_intervals_NI = [bite_intervals_NI diff(bite_timestamps_temp)];
                        end
                    end
                else
                    trials_NI(end+1).bite_timestamps = bite_timestamps;
                    trials_NI(end).bite_amplitudes = bite_amplitudes;
                    if ~isempty(bite_timestamps)
                        firstbite_time_NI(end+1) = bite_timestamps(1);
                    end
                    totalbites_NI(end+1) = numel(bite_timestamps);
                    if app.BiteDeviceButton.Value
                        bite_duration_NI(end+1) = 0;
                        for j = 1:size(Bite_events(i).bite_window, 1)
                            bite_timestamps_temp = bite_timestamps(bite_timestamps > Bite_events(i).bite_window(j, 1) & bite_timestamps < Bite_events(i).bite_window(j, 2));
                            bite_duration_NI(end) = bite_duration_NI(end)+range(bite_timestamps_temp);
                            bite_intervals_NI = [bite_intervals_NI diff(bite_timestamps_temp)];
                        end
                    end
                end
                if ~app.BiteDeviceButton.Value
                    bite_intervals_NI = [bite_intervals_NI diff(bite_timestamps)];
                end
                bite_amplitudes_NI = [bite_amplitudes_NI bite_amplitudes];
            else
                if isempty(trials_I)
                    trials_I(1).bite_timestamps = bite_timestamps;
                    trials_I(1).bite_amplitudes = bite_amplitudes;
                    trials_I(1).laser_timestamps = Bite_events(i).laser_timestamps;
                    if ~isempty(bite_timestamps)
                        firstbite_time_I(1) = bite_timestamps(1);
                    end
                    totalbites_I(1) = numel(find(bite_timestamps >= Bite_events(i).laser_timestamps(1) &...
                        bite_timestamps <= Bite_events(i).laser_timestamps(2)));
                    if app.BiteDeviceButton.Value
                        bite_duration_I(1) = 0;
                        for j = 1:size(Bite_events(i).bite_window, 1)
                            bite_timestamps_temp = bite_timestamps(bite_timestamps > Bite_events(i).bite_window(j, 1) & bite_timestamps < Bite_events(i).bite_window(j, 2));
                            bite_duration_I(1) = bite_duration_I(1)+range(bite_timestamps_temp);
                            bite_intervals_I = [bite_intervals_I diff(bite_timestamps_temp)];
                        end
                    end
                else
                    trials_I(end+1).bite_timestamps = bite_timestamps;
                    trials_I(end).bite_amplitudes = bite_amplitudes;
                    trials_I(end).laser_timestamps = Bite_events(i).laser_timestamps;
                    if ~isempty(bite_timestamps)
                        firstbite_time_I(end+1) = bite_timestamps(1);
                    end
                    totalbites_I(end+1) = numel(find(bite_timestamps >= Bite_events(i).laser_timestamps(1) &...
                        bite_timestamps <= Bite_events(i).laser_timestamps(2)));
                    if app.BiteDeviceButton.Value
                        bite_duration_I(end+1) = 0;
                        for j = 1:size(Bite_events(i).bite_window, 1)
                            bite_timestamps_temp = bite_timestamps(bite_timestamps > Bite_events(i).bite_window(j, 1) & bite_timestamps < Bite_events(i).bite_window(j, 2));
                            bite_duration_I(end) = bite_duration_I(end)+range(bite_timestamps_temp);
                            bite_intervals_I = [bite_intervals_I diff(bite_timestamps_temp)];
                        end
                    end
                end
                if ~app.BiteDeviceButton.Value
                    bite_intervals_I = [bite_intervals_I diff(bite_timestamps(bite_timestamps >= Bite_events(i).laser_timestamps(1) &...
                        bite_timestamps <= Bite_events(i).laser_timestamps(2)))];
                end
                bite_amplitudes_I = [bite_amplitudes_I bite_amplitudes(bite_timestamps >= Bite_events(i).laser_timestamps(1) &...
                    bite_timestamps <= Bite_events(i).laser_timestamps(2))];
            end
        elseif app.sOnTwiceButton.Value
            if isempty(Bite_events(i).laser_timestamps)
                if isempty(trials_NI)
                    trials_NI(1).bite_timestamps = bite_timestamps;
                    trials_NI(1).bite_amplitudes = bite_amplitudes;
                    if ~isempty(bite_timestamps)
                        firstafter2lastbefore_NI1(1) = bite_timestamps(find(bite_timestamps >= 4, 1))-...
                            bite_timestamps(find(bite_timestamps < 4, 1, 'last'));
                        firstafter2lastbefore_NI2(1) = bite_timestamps(find(bite_timestamps >= 13, 1))-...
                            bite_timestamps(find(bite_timestamps < 13, 1, 'last'));
                    end
                    totalbites_NI1(1) = numel(find(bite_timestamps >= 4 & bite_timestamps <= 8));
                    totalbites_NI2(1) = numel(find(bite_timestamps >= 13 & bite_timestamps <= 17));
                else
                    trials_NI(end+1).bite_timestamps = bite_timestamps;
                    trials_NI(end).bite_amplitudes = bite_amplitudes;
                    if ~isempty(bite_timestamps)
                        firstafter2lastbefore_NI1(end+1) = bite_timestamps(find(bite_timestamps >= 4, 1))-...
                            bite_timestamps(find(bite_timestamps < 4, 1, 'last'));
                        firstafter2lastbefore_NI2(end+1) = bite_timestamps(find(bite_timestamps >= 13, 1))-...
                            bite_timestamps(find(bite_timestamps < 13, 1, 'last'));
                    end
                    totalbites_NI1(end+1) = numel(find(bite_timestamps >= 4 & bite_timestamps <= 8));
                    totalbites_NI2(end+1) = numel(find(bite_timestamps >= 13 & bite_timestamps <= 17));
                end
                bite_intervals_NI1 = [bite_intervals_NI1 diff(bite_timestamps(bite_timestamps >= 4 & bite_timestamps <= 8))];
                bite_amplitudes_NI1 = [bite_amplitudes_NI1 bite_amplitudes(bite_timestamps >= 4 & bite_timestamps <= 8)];
                bite_intervals_NI2 = [bite_intervals_NI2 diff(bite_timestamps(bite_timestamps >= 13 & bite_timestamps <= 17))];
                bite_amplitudes_NI2 = [bite_amplitudes_NI2 bite_amplitudes(bite_timestamps >= 13 & bite_timestamps <= 17)];
            else
                if isempty(trials_I)
                    trials_I(1).bite_timestamps = bite_timestamps;
                    trials_I(1).bite_amplitudes = bite_amplitudes;
                    trials_I(1).laser_timestamps = Bite_events(i).laser_timestamps;
                    if ~isempty(bite_timestamps)
                        firstafter2lastbefore_I1(1) = bite_timestamps(find(bite_timestamps >= Bite_events(i).laser_timestamps(1, 1), 1))-...
                            bite_timestamps(find(bite_timestamps < Bite_events(i).laser_timestamps(1, 1), 1, 'last'));
                        firstafter2lastbefore_I2(1) = bite_timestamps(find(bite_timestamps >= Bite_events(i).laser_timestamps(1, 2), 1))-...
                            bite_timestamps(find(bite_timestamps < Bite_events(i).laser_timestamps(1, 2), 1, 'last'));
                    end
                    totalbites_I1(1) = numel(find(bite_timestamps >= Bite_events(i).laser_timestamps(1, 1) &...
                        bite_timestamps <= Bite_events(i).laser_timestamps(2, 1)));
                    totalbites_I2(1) = numel(find(bite_timestamps >= Bite_events(i).laser_timestamps(1, 2) &...
                        bite_timestamps <= Bite_events(i).laser_timestamps(2, 2)));
                else
                    trials_I(end+1).bite_timestamps = bite_timestamps;
                    trials_I(end).bite_amplitudes = bite_amplitudes;
                    trials_I(end).laser_timestamps = Bite_events(i).laser_timestamps;
                    if ~isempty(bite_timestamps)
                        firstafter2lastbefore_I1(end+1) = bite_timestamps(find(bite_timestamps >= Bite_events(i).laser_timestamps(1, 1), 1))-...
                            bite_timestamps(find(bite_timestamps < Bite_events(i).laser_timestamps(1, 1), 1, 'last'));
                        firstafter2lastbefore_I2(end+1) = bite_timestamps(find(bite_timestamps >= Bite_events(i).laser_timestamps(1, 2), 1))-...
                            bite_timestamps(find(bite_timestamps < Bite_events(i).laser_timestamps(1, 2), 1, 'last'));
                    end
                    totalbites_I1(end+1) = numel(find(bite_timestamps >= Bite_events(i).laser_timestamps(1, 1) &...
                        bite_timestamps <= Bite_events(i).laser_timestamps(2, 1)));
                    totalbites_I2(end+1) = numel(find(bite_timestamps >= Bite_events(i).laser_timestamps(1, 2) &...
                        bite_timestamps <= Bite_events(i).laser_timestamps(2, 2)));
                end
                bite_intervals_I1 = [bite_intervals_I1 diff(bite_timestamps(bite_timestamps >= Bite_events(i).laser_timestamps(1, 1) &...
                    bite_timestamps <= Bite_events(i).laser_timestamps(2, 1)))];
                bite_amplitudes_I1 = [bite_amplitudes_I1 bite_amplitudes(bite_timestamps >= Bite_events(i).laser_timestamps(1, 1) &...
                    bite_timestamps <= Bite_events(i).laser_timestamps(2, 1))];
                bite_intervals_I2 = [bite_intervals_I2 diff(bite_timestamps(bite_timestamps >= Bite_events(i).laser_timestamps(1, 2) &...
                    bite_timestamps <= Bite_events(i).laser_timestamps(2, 2)))];
                bite_amplitudes_I2 = [bite_amplitudes_I2 bite_amplitudes(bite_timestamps >= Bite_events(i).laser_timestamps(1, 2) &...
                    bite_timestamps <= Bite_events(i).laser_timestamps(2, 2))];
            end
        end
    end
end

if app.NoLightButton.Value
    bite_crosscorrelogram(all_trials);
elseif app.WholeTrialOnButton.Value || app.BiteDeviceButton.Value
    bite_crosscorrelogram(trials_NI, trials_I);
else
    bite_crosscorrelogram(trials_NI);
end

time_range = [app.timerangesEditField.Value app.toEditField.Value];
bite_edges = 0:bite_bin:max_timestamp;
if app.sDelay4sOnButton.Value
    bite_rate_NI = zeros(numel(trials_NI), numel(bite_edges)-1);
    bite_rate_I = zeros(numel(trials_I), numel(bite_edges)-1);
    % raster plot
    figure;
    hold on;
    for i = 1:numel(trials_NI)
        bite_rate_NI(i, :) = histcounts(trials_NI(i).bite_timestamps, bite_edges)/bite_bin;
        for j = 1:numel(trials_NI(i).bite_timestamps)
            line([trials_NI(i).bite_timestamps(j); trials_NI(i).bite_timestamps(j)],...
                [i-0.5-trials_NI(i).bite_amplitudes(j)/2; i-0.5+trials_NI(i).bite_amplitudes(j)/2],...
                'Color', 'k', 'LineStyle', '-', 'LineWidth', 1);
        end
    end
    line([4 8; 4 8], [0 0; numel(trials_NI) numel(trials_NI)], 'Color', 'r', 'LineStyle', '--', 'LineWidth', 1);
    for i = 1:numel(trials_I)
        bite_rate_I(i, :) = histcounts(trials_I(i).bite_timestamps, bite_edges)/bite_bin;
        patch([trials_I(i).laser_timestamps(1) trials_I(i).laser_timestamps(1) trials_I(i).laser_timestamps(2) trials_I(i).laser_timestamps(2)],...
            [numel(trials_NI)+i-1 numel(trials_NI)+i numel(trials_NI)+i numel(trials_NI)+i-1], [0 1 0], 'FaceAlpha', 1, 'EdgeColor', 'none');
        for j = 1:numel(trials_I(i).bite_timestamps)
            line([trials_I(i).bite_timestamps(j); trials_I(i).bite_timestamps(j)],...
                [numel(trials_NI)+i-0.5-trials_I(i).bite_amplitudes(j)/2; numel(trials_NI)+i-0.5+trials_I(i).bite_amplitudes(j)/2],...
                'Color', 'k', 'LineStyle', '-', 'LineWidth', 1);
        end
    end
    if ~isinf(time_range(2))
        xlim(time_range);
    elseif isinf(time_range(2))
        xl = xlim;
        xl(1) = time_range(1);
        xlim(xl);
    end
    set(gca, 'YTick', 0.5:1:numel(trials_NI)+numel(trials_I)-0.5, 'YTickLabel', num2str((1:1:numel(trials_NI)+numel(trials_I))'), 'yLim', [0 numel(trials_I)+numel(trials_NI)],...
        'TickLength', [0 0], 'FontSize', 12, 'YDir', 'reverse');
    xlabel('Time (s)');
    ylabel('Trials');
    
    figure;
    plot_tj_multitrial(time_range, bite_edges(2:end)', bite_rate_NI', bite_edges(2:end)', bite_rate_I',...
        get_laser_timestamp_sound(trials_I, 'start', 1, 'mean'), get_laser_timestamp_sound(trials_I, 'stop', 1, 'mean'), 'Bite (bite/s)', '');
    
    % bar plots
    figure;
    subplot(1, 4, 1);
    hold on;
    bar(1, mean(firstafter2lastbefore_NI), 0.8, 'FaceColor', [0.5 0.5 0.5]);
    plot(1+rand(1, numel(firstafter2lastbefore_NI))*0.8-0.4, firstafter2lastbefore_NI, 'o',...
        'MarkerFaceColor', [0 0 0], 'MarkerEdgeColor', 'none');
    bar(2, mean(firstafter2lastbefore_I), 0.8, 'FaceColor', [0 1 0]);
    plot(2+rand(1, numel(firstafter2lastbefore_I))*0.8-0.4, firstafter2lastbefore_I, 'o',...
        'MarkerFaceColor', [0 0.75 0], 'MarkerEdgeColor', 'none');
    errorbar([1 2], [mean(firstafter2lastbefore_NI) mean(firstafter2lastbefore_I)],...
        [std(firstafter2lastbefore_NI)/sqrt(numel(firstafter2lastbefore_NI)) std(firstafter2lastbefore_I)/sqrt(numel(firstafter2lastbefore_I))],...
        'k', 'LineStyle', 'none');
    set(gca, 'XTick', [1 2], 'XTickLabel', {'Control', 'Inhibition'}, 'XLim', [0.4 2.6], 'TickLength', [0 0], 'FontSize', 12);
    ylabel('Interval (s)');
    fprintf('\n');
    p = ranksum(firstafter2lastbefore_NI, firstafter2lastbefore_I, 'tail', 'both');
    disp(['firstafter2lastbefore (both): p = ' num2str(p)]);
    p = ranksum(firstafter2lastbefore_NI, firstafter2lastbefore_I, 'tail', 'left');
    disp(['firstafter2lastbefore (left): p = ' num2str(p)]);
    
    subplot(1, 4, 2);
    hold on;
    bar(1, mean(totalbites_NI), 0.8, 'FaceColor', [0.5 0.5 0.5]);
    plot(1+rand(1, numel(totalbites_NI))*0.8-0.4, totalbites_NI, 'o', 'MarkerFaceColor', [0 0 0], 'MarkerEdgeColor', 'none');
    bar(2, mean(totalbites_I), 0.8, 'FaceColor', [0 1 0]);
    plot(2+rand(1, numel(totalbites_I))*0.8-0.4, totalbites_I, 'o', 'MarkerFaceColor', [0 0.75 0], 'MarkerEdgeColor', 'none');
    errorbar([1 2], [mean(totalbites_NI) mean(totalbites_I)],...
        [std(totalbites_NI)/sqrt(numel(totalbites_NI)) std(totalbites_I)/sqrt(numel(totalbites_I))], 'k', 'LineStyle', 'none');
    set(gca, 'XTick', [1 2], 'XTickLabel', {'Control', 'Inhibition'}, 'XLim', [0.4 2.6], 'TickLength', [0 0], 'FontSize', 12);
    ylabel('Bites');
    fprintf('\n');
    p = ranksum(totalbites_NI, totalbites_I, 'tail', 'both');
    disp(['totalbites (both): p = ' num2str(p)]);
    p = ranksum(totalbites_NI, totalbites_I, 'tail', 'right');
    disp(['totalbites (right): p = ' num2str(p)]);
    
    try
        subplot(1, 4, 3);
        hold on;
        bar(1, mean(bite_intervals_NI), 0.8, 'FaceColor', [0.5 0.5 0.5]);
        plot(1+rand(1, numel(bite_intervals_NI))*0.8-0.4, bite_intervals_NI, 'o', 'MarkerFaceColor', [0 0 0], 'MarkerEdgeColor', 'none');
        bar(2, mean(bite_intervals_I), 0.8, 'FaceColor', [0 1 0]);
        plot(2+rand(1, numel(bite_intervals_I))*0.8-0.4, bite_intervals_I, 'o', 'MarkerFaceColor', [0 0.75 0], 'MarkerEdgeColor', 'none');
        errorbar([1 2], [mean(bite_intervals_NI) mean(bite_intervals_I)],...
            [std(bite_intervals_NI)/sqrt(numel(bite_intervals_NI)) std(bite_intervals_I)/sqrt(numel(bite_intervals_I))], 'k', 'LineStyle', 'none');
        set(gca, 'XTick', [1 2], 'XTickLabel', {'Control', 'Inhibition'}, 'XLim', [0.4 2.6], 'TickLength', [0 0], 'FontSize', 12);
        ylabel('Bite Interval (s)');
        fprintf('\n');
        p = ranksum(bite_intervals_NI, bite_intervals_I, 'tail', 'both');
        disp(['bite interval (both): p = ' num2str(p)]);
        p = ranksum(bite_intervals_NI, bite_intervals_I, 'tail', 'left');
        disp(['bite interval (left): p = ' num2str(p)]);
    end
    try
        subplot(1, 4, 4);
        hold on;
        bar(1, mean(bite_amplitudes_NI), 0.8, 'FaceColor', [0.5 0.5 0.5]);
        plot(1+rand(1, numel(bite_amplitudes_NI))*0.8-0.4, bite_amplitudes_NI, 'o', 'MarkerFaceColor', [0 0 0], 'MarkerEdgeColor', 'none');
        bar(2, mean(bite_amplitudes_I), 0.8, 'FaceColor', [0 1 0]);
        plot(2+rand(1, numel(bite_amplitudes_I))*0.8-0.4, bite_amplitudes_I, 'o', 'MarkerFaceColor', [0 0.75 0], 'MarkerEdgeColor', 'none');
        errorbar([1 2], [mean(bite_amplitudes_NI) mean(bite_amplitudes_I)],...
            [std(bite_amplitudes_NI)/sqrt(numel(bite_amplitudes_NI)) std(bite_amplitudes_I)/sqrt(numel(bite_amplitudes_I))],...
            'k', 'LineStyle', 'none');
        set(gca, 'XTick', [1 2], 'XTickLabel', {'Control', 'Inhibition'}, 'XLim', [0.4 2.6], 'TickLength', [0 0], 'FontSize', 12);
        ylabel('Normalized Sound Amplitude');
        fprintf('\n');
        p = ranksum(bite_amplitudes_NI, bite_amplitudes_I, 'tail', 'both');
        disp(['bite amplitude (both): p = ' num2str(p)]);
        p = ranksum(bite_amplitudes_NI, bite_amplitudes_I, 'tail', 'right');
        disp(['bite amplitude (right): p = ' num2str(p)]);
    end
elseif app.sDelay4sOnButton_2.Value
    bite_rate_NI = zeros(numel(trials_NI), numel(bite_edges)-1);
    bite_rate_I = zeros(numel(trials_I), numel(bite_edges)-1);
    % raster plot
    figure;
    hold on;
    for i = 1:numel(trials_NI)
        bite_rate_NI(i, :) = histcounts(trials_NI(i).bite_timestamps, bite_edges)/bite_bin;
        for j = 1:numel(trials_NI(i).bite_timestamps)
            line([trials_NI(i).bite_timestamps(j); trials_NI(i).bite_timestamps(j)],...
                [i-0.5-trials_NI(i).bite_amplitudes(j)/2; i-0.5+trials_NI(i).bite_amplitudes(j)/2],...
                'Color', 'k', 'LineStyle', '-', 'LineWidth', 1);
        end
    end
    line([4; 4], [0; numel(trials_NI)], 'Color', 'r', 'LineStyle', '--', 'LineWidth', 1);
    for i = 1:numel(trials_I)
        bite_rate_I(i, :) = histcounts(trials_I(i).bite_timestamps, bite_edges)/bite_bin;
        patch([trials_I(i).laser_timestamps(1) trials_I(i).laser_timestamps(1) trials_I(i).laser_timestamps(2) trials_I(i).laser_timestamps(2)],...
            [numel(trials_NI)+i-1 numel(trials_NI)+i numel(trials_NI)+i numel(trials_NI)+i-1], [0 1 0], 'FaceAlpha', 1, 'EdgeColor', 'none');
        for j = 1:numel(trials_I(i).bite_timestamps)
            line([trials_I(i).bite_timestamps(j); trials_I(i).bite_timestamps(j)],...
                [numel(trials_NI)+i-0.5-trials_I(i).bite_amplitudes(j)/2; numel(trials_NI)+i-0.5+trials_I(i).bite_amplitudes(j)/2],...
                'Color', 'k', 'LineStyle', '-', 'LineWidth', 1);
        end
    end
    if ~isinf(time_range(2))
        xlim(time_range);
    elseif isinf(time_range(2))
        xl = xlim;
        xl(1) = time_range(1);
        xlim(xl);
    end
    set(gca, 'YTick', 0.5:1:numel(trials_NI)+numel(trials_I)-0.5, 'YTickLabel', num2str((1:1:numel(trials_NI)+numel(trials_I))'), 'yLim', [0 numel(trials_I)+numel(trials_NI)],...
        'TickLength', [0 0], 'FontSize', 12, 'YDir', 'reverse');
    xlabel('Time (s)');
    ylabel('Trials');
    
    figure;
    plot_tj_multitrial(time_range, bite_edges(2:end)', bite_rate_NI', bite_edges(2:end)', bite_rate_I',...
        get_laser_timestamp_sound(trials_I, 'start', 1, 'mean'), get_laser_timestamp_sound(trials_I, 'stop', 1, 'mean'), 'Bite (bite/s)', '');
    
    % bar plots
    figure;
    subplot(1, 4, 1);
    hold on;
    bar(1, mean(firstbite_time_NI), 0.8, 'FaceColor', [0.5 0.5 0.5]);
    plot(1+rand(1, numel(firstbite_time_NI))*0.8-0.4, firstbite_time_NI, 'o', 'MarkerFaceColor', [0 0 0], 'MarkerEdgeColor', 'none');
    bar(2, mean(firstbite_time_I), 0.8, 'FaceColor', [0 1 0]);
    plot(2+rand(1, numel(firstbite_time_I))*0.8-0.4, firstbite_time_I, 'o', 'MarkerFaceColor', [0 0.75 0], 'MarkerEdgeColor', 'none');
    errorbar([1 2], [mean(firstbite_time_NI) mean(firstbite_time_I)],...
        [std(firstbite_time_NI)/sqrt(numel(firstbite_time_NI)) std(firstbite_time_I)/sqrt(numel(firstbite_time_I))],...
        'k', 'LineStyle', 'none');
    set(gca, 'XTick', [1 2], 'XTickLabel', {'Control', 'Inhibition'}, 'XLim', [0.4 2.6], 'TickLength', [0 0], 'FontSize', 12);
    ylabel('First Bite (s)');
    fprintf('\n');
    p = ranksum(firstbite_time_NI, firstbite_time_I, 'tail', 'both');
    disp(['First Bite (both): p = ' num2str(p)]);
    p = ranksum(firstbite_time_NI, firstbite_time_I, 'tail', 'left');
    disp(['First Bite (left): p = ' num2str(p)]);
    
    subplot(1, 4, 2);
    hold on;
    bar(1, mean(totalbites_NI), 0.8, 'FaceColor', [0.5 0.5 0.5]);
    plot(1+rand(1, numel(totalbites_NI))*0.8-0.4, totalbites_NI, 'o', 'MarkerFaceColor', [0 0 0], 'MarkerEdgeColor', 'none');
    bar(2, mean(totalbites_I), 0.8, 'FaceColor', [0 1 0]);
    plot(2+rand(1, numel(totalbites_I))*0.8-0.4, totalbites_I, 'o', 'MarkerFaceColor', [0 0.75 0], 'MarkerEdgeColor', 'none');
    errorbar([1 2], [mean(totalbites_NI) mean(totalbites_I)],...
        [std(totalbites_NI)/sqrt(numel(totalbites_NI)) std(totalbites_I)/sqrt(numel(totalbites_I))], 'k', 'LineStyle', 'none');
    set(gca, 'XTick', [1 2], 'XTickLabel', {'Control', 'Inhibition'}, 'XLim', [0.4 2.6], 'TickLength', [0 0], 'FontSize', 12);
    ylabel('Bites');
    fprintf('\n');
    p = ranksum(totalbites_NI, totalbites_I, 'tail', 'both');
    disp(['totalbites (both): p = ' num2str(p)]);
    p = ranksum(totalbites_NI, totalbites_I, 'tail', 'right');
    disp(['totalbites (right): p = ' num2str(p)]);
    
    try
        subplot(1, 4, 3);
        hold on;
        bar(1, mean(bite_intervals_NI), 0.8, 'FaceColor', [0.5 0.5 0.5]);
        plot(1+rand(1, numel(bite_intervals_NI))*0.8-0.4, bite_intervals_NI, 'o', 'MarkerFaceColor', [0 0 0], 'MarkerEdgeColor', 'none');
        bar(2, mean(bite_intervals_I), 0.8, 'FaceColor', [0 1 0]);
        plot(2+rand(1, numel(bite_intervals_I))*0.8-0.4, bite_intervals_I, 'o', 'MarkerFaceColor', [0 0.75 0], 'MarkerEdgeColor', 'none');
        errorbar([1 2], [mean(bite_intervals_NI) mean(bite_intervals_I)],...
            [std(bite_intervals_NI)/sqrt(numel(bite_intervals_NI)) std(bite_intervals_I)/sqrt(numel(bite_intervals_I))], 'k', 'LineStyle', 'none');
        set(gca, 'XTick', [1 2], 'XTickLabel', {'Control', 'Inhibition'}, 'XLim', [0.4 2.6], 'TickLength', [0 0], 'FontSize', 12);
        ylabel('Bite Interval (s)');
        fprintf('\n');
        p = ranksum(bite_intervals_NI, bite_intervals_I, 'tail', 'both');
        disp(['bite interval (both): p = ' num2str(p)]);
        p = ranksum(bite_intervals_NI, bite_intervals_I, 'tail', 'left');
        disp(['bite interval (left): p = ' num2str(p)]);
    end
    
    try
        subplot(1, 4, 4);
        hold on;
        bar(1, mean(bite_amplitudes_NI), 0.8, 'FaceColor', [0.5 0.5 0.5]);
        plot(1+rand(1, numel(bite_amplitudes_NI))*0.8-0.4, bite_amplitudes_NI, 'o', 'MarkerFaceColor', [0 0 0], 'MarkerEdgeColor', 'none');
        bar(2, mean(bite_amplitudes_I), 0.8, 'FaceColor', [0 1 0]);
        plot(2+rand(1, numel(bite_amplitudes_I))*0.8-0.4, bite_amplitudes_I, 'o', 'MarkerFaceColor', [0 0.75 0], 'MarkerEdgeColor', 'none');
        errorbar([1 2], [mean(bite_amplitudes_NI) mean(bite_amplitudes_I)],...
            [std(bite_amplitudes_NI)/sqrt(numel(bite_amplitudes_NI)) std(bite_amplitudes_I)/sqrt(numel(bite_amplitudes_I))],...
            'k', 'LineStyle', 'none');
        set(gca, 'XTick', [1 2], 'XTickLabel', {'Control', 'Inhibition'}, 'XLim', [0.4 2.6], 'TickLength', [0 0], 'FontSize', 12);
        ylabel('Normalized Sound Amplitude');
        fprintf('\n');
        p = ranksum(bite_amplitudes_NI, bite_amplitudes_I, 'tail', 'both');
        disp(['bite amplitude (both): p = ' num2str(p)]);
        p = ranksum(bite_amplitudes_NI, bite_amplitudes_I, 'tail', 'right');
        disp(['bite amplitude (right): p = ' num2str(p)]);
    end
elseif app.NoLightButton.Value
    bite_rate = zeros(numel(all_trials), numel(bite_edges)-1);
    % raster plot
    figure;
    hold on;
    for i = 1:numel(all_trials)
        bite_rate(i, :) = histcounts(all_trials(i).bite_timestamps, bite_edges)/bite_bin;
        for j = 1:numel(all_trials(i).bite_timestamps)
            line([all_trials(i).bite_timestamps(j); all_trials(i).bite_timestamps(j)],...
                [i-0.5-all_trials(i).bite_amplitudes(j)/2; i-0.5+all_trials(i).bite_amplitudes(j)/2],...
                'Color', 'k', 'LineStyle', '-', 'LineWidth', 1);
        end
    end
    if ~isinf(time_range(2))
        xlim(time_range);
    elseif isinf(time_range(2))
        xl = xlim;
        xl(1) = time_range(1);
        xlim(xl);
    end
    set(gca, 'YTick', 0.5:1:numel(all_trials)-0.5, 'YTickLabel', num2str((1:1:numel(all_trials))'), 'yLim', [0 numel(all_trials)], 'TickLength', [0 0], 'FontSize', 12, 'YDir', 'reverse');
    xlabel('Time (s)');
    ylabel('Trials');
    
    figure;
    plot_tj_multitrial(time_range, bite_edges(2:end)', bite_rate', [], [],...
        [], [], 'Bite (bite/s)', '');
    
    % bar plots
    figure;
    subplot(1, 3, 1);
    hold on;
    bar(1, mean(firstbite_time), 0.8, 'FaceColor', [0.5 0.5 0.5]);
    plot(1+rand(1, numel(firstbite_time))*0.8-0.4, firstbite_time, 'o', 'MarkerFaceColor', [0 0 0], 'MarkerEdgeColor', 'none');
    errorbar(1, mean(firstbite_time), std(firstbite_time)/sqrt(numel(firstbite_time)), 'k', 'LineStyle', 'none');
    set(gca, 'XTick', 1, 'XTickLabel', {''}, 'XLim', [0.4 1.6], 'TickLength', [0 0], 'FontSize', 12);
    ylabel('First Bite (s)');
    
    subplot(1, 3, 2);
    hold on;
    bar(1, mean(totalbites), 0.8, 'FaceColor', [0.5 0.5 0.5]);
    plot(1+rand(1, numel(totalbites))*0.8-0.4, totalbites, 'o', 'MarkerFaceColor', [0 0 0], 'MarkerEdgeColor', 'none');
    errorbar(1, mean(totalbites), std(totalbites)/sqrt(numel(totalbites)), 'k', 'LineStyle', 'none');
    set(gca, 'XTick', 1, 'XTickLabel', {''}, 'XLim', [0.4 1.6], 'TickLength', [0 0], 'FontSize', 12);
    ylabel('Bites');
    
    subplot(1, 3, 3);
    hold on;
    bar(1, mean(bite_intervals), 0.8, 'FaceColor', [0.5 0.5 0.5]);
    plot(1+rand(1, numel(bite_intervals))*0.8-0.4, bite_intervals, 'o', 'MarkerFaceColor', [0 0 0], 'MarkerEdgeColor', 'none');
    errorbar(1, mean(bite_intervals), std(bite_intervals)/sqrt(numel(bite_intervals)), 'k', 'LineStyle', 'none');
    set(gca, 'XTick', 1, 'XTickLabel', {''}, 'XLim', [0.4 1.6], 'TickLength', [0 0], 'FontSize', 12);
    ylabel('Bite Interval (s)');
elseif app.WholeTrialOnButton.Value || app.BiteDeviceButton.Value
    bite_rate_NI = zeros(numel(trials_NI), numel(bite_edges)-1);
    bite_rate_I = zeros(numel(trials_I), numel(bite_edges)-1);
    % raster plot
    figure;
    hold on;
    for i = 1:numel(trials_NI)
        bite_rate_NI(i, :) = histcounts(trials_NI(i).bite_timestamps, bite_edges)/bite_bin;
        for j = 1:numel(trials_NI(i).bite_timestamps)
            line([trials_NI(i).bite_timestamps(j); trials_NI(i).bite_timestamps(j)],...
                [i-0.5-trials_NI(i).bite_amplitudes(j)/2; i-0.5+trials_NI(i).bite_amplitudes(j)/2],...
                'Color', 'k', 'LineStyle', '-', 'LineWidth', 1);
        end
    end
    for i = 1:numel(trials_I)
        bite_rate_I(i, :) = histcounts(trials_I(i).bite_timestamps, bite_edges)/bite_bin;
        patch([trials_I(i).laser_timestamps(1) trials_I(i).laser_timestamps(1) trials_I(i).laser_timestamps(2) trials_I(i).laser_timestamps(2)],...
            [numel(trials_NI)+i-1 numel(trials_NI)+i numel(trials_NI)+i numel(trials_NI)+i-1], [0 1 0], 'FaceAlpha', 1, 'EdgeColor', 'none');
        for j = 1:numel(trials_I(i).bite_timestamps)
            line([trials_I(i).bite_timestamps(j); trials_I(i).bite_timestamps(j)],...
                [numel(trials_NI)+i-0.5-trials_I(i).bite_amplitudes(j)/2; numel(trials_NI)+i-0.5+trials_I(i).bite_amplitudes(j)/2],...
                'Color', 'k', 'LineStyle', '-', 'LineWidth', 1);
        end
    end
    if ~isinf(time_range(2))
        xlim(time_range);
    elseif isinf(time_range(2))
        xl = xlim;
        xl(1) = time_range(1);
        xlim(xl);
    end
    set(gca, 'YTick', 0.5:1:numel(trials_NI)+numel(trials_I)-0.5, 'YTickLabel', num2str((1:1:numel(trials_NI)+numel(trials_I))'), 'yLim', [0 numel(trials_I)+numel(trials_NI)],...
        'TickLength', [0 0], 'FontSize', 12, 'YDir', 'reverse');
    xlabel('Time (s)');
    ylabel('Trials');
    
    figure;
    plot_tj_multitrial(time_range, bite_edges(2:end)', bite_rate_NI', bite_edges(2:end)', bite_rate_I',...
        [], [], 'Bite (bite/s)', '');
    
    % bar plots
    figure;
    subplot(1, 4, 1);
    hold on;
    bar(1, mean(firstbite_time_NI), 0.8, 'FaceColor', [0.5 0.5 0.5]);
    plot(1+rand(1, numel(firstbite_time_NI))*0.8-0.4, firstbite_time_NI, 'o', 'MarkerFaceColor', [0 0 0], 'MarkerEdgeColor', 'none');
    bar(2, mean(firstbite_time_I), 0.8, 'FaceColor', [0 1 0]);
    plot(2+rand(1, numel(firstbite_time_I))*0.8-0.4, firstbite_time_I, 'o', 'MarkerFaceColor', [0 0.75 0], 'MarkerEdgeColor', 'none');
    errorbar([1 2], [mean(firstbite_time_NI) mean(firstbite_time_I)],...
        [std(firstbite_time_NI)/sqrt(numel(firstbite_time_NI)) std(firstbite_time_I)/sqrt(numel(firstbite_time_I))],...
        'k', 'LineStyle', 'none');
    set(gca, 'XTick', [1 2], 'XTickLabel', {'Control', 'Inhibition'}, 'XLim', [0.4 2.6], 'TickLength', [0 0], 'FontSize', 12);
    ylabel('First Bite (s)');
    fprintf('\n');
    p = ranksum(firstbite_time_NI, firstbite_time_I, 'tail', 'both');
    disp(['First Bite (both): p = ' num2str(p)]);
    p = ranksum(firstbite_time_NI, firstbite_time_I, 'tail', 'left');
    disp(['First Bite (left): p = ' num2str(p)]);
    
    subplot(1, 4, 2);
    hold on;
    bar(1, mean(totalbites_NI), 0.8, 'FaceColor', [0.5 0.5 0.5]);
    plot(1+rand(1, numel(totalbites_NI))*0.8-0.4, totalbites_NI, 'o', 'MarkerFaceColor', [0 0 0], 'MarkerEdgeColor', 'none');
    bar(2, mean(totalbites_I), 0.8, 'FaceColor', [0 1 0]);
    plot(2+rand(1, numel(totalbites_I))*0.8-0.4, totalbites_I, 'o', 'MarkerFaceColor', [0 0.75 0], 'MarkerEdgeColor', 'none');
    errorbar([1 2], [mean(totalbites_NI) mean(totalbites_I)],...
        [std(totalbites_NI)/sqrt(numel(totalbites_NI)) std(totalbites_I)/sqrt(numel(totalbites_I))], 'k', 'LineStyle', 'none');
    set(gca, 'XTick', [1 2], 'XTickLabel', {'Control', 'Inhibition'}, 'XLim', [0.4 2.6], 'TickLength', [0 0], 'FontSize', 12);
    ylabel('Bites');
    fprintf('\n');
    p = ranksum(totalbites_NI, totalbites_I, 'tail', 'both');
    disp(['totalbites (both): p = ' num2str(p)]);
    p = ranksum(totalbites_NI, totalbites_I, 'tail', 'right');
    disp(['totalbites (right): p = ' num2str(p)]);
    
    subplot(1, 4, 3);
    hold on;
    bar(1, mean(bite_intervals_NI), 0.8, 'FaceColor', [0.5 0.5 0.5]);
    plot(1+rand(1, numel(bite_intervals_NI))*0.8-0.4, bite_intervals_NI, 'o', 'MarkerFaceColor', [0 0 0], 'MarkerEdgeColor', 'none');
    bar(2, mean(bite_intervals_I), 0.8, 'FaceColor', [0 1 0]);
    plot(2+rand(1, numel(bite_intervals_I))*0.8-0.4, bite_intervals_I, 'o', 'MarkerFaceColor', [0 0.75 0], 'MarkerEdgeColor', 'none');
    errorbar([1 2], [mean(bite_intervals_NI) mean(bite_intervals_I)],...
        [std(bite_intervals_NI)/sqrt(numel(bite_intervals_NI)) std(bite_intervals_I)/sqrt(numel(bite_intervals_I))], 'k', 'LineStyle', 'none');
    set(gca, 'XTick', [1 2], 'XTickLabel', {'Control', 'Inhibition'}, 'XLim', [0.4 2.6], 'TickLength', [0 0], 'FontSize', 12);
    ylabel('Bite Interval (s)');
    try
        fprintf('\n');
        p = ranksum(bite_intervals_NI, bite_intervals_I, 'tail', 'both');
        disp(['bite interval (both): p = ' num2str(p)]);
        p = ranksum(bite_intervals_NI, bite_intervals_I, 'tail', 'left');
        disp(['bite interval (left): p = ' num2str(p)]);
    end
    
    subplot(1, 4, 4);
    hold on;
    bar(1, mean(bite_amplitudes_NI), 0.8, 'FaceColor', [0.5 0.5 0.5]);
    plot(1+rand(1, numel(bite_amplitudes_NI))*0.8-0.4, bite_amplitudes_NI, 'o', 'MarkerFaceColor', [0 0 0], 'MarkerEdgeColor', 'none');
    bar(2, mean(bite_amplitudes_I), 0.8, 'FaceColor', [0 1 0]);
    plot(2+rand(1, numel(bite_amplitudes_I))*0.8-0.4, bite_amplitudes_I, 'o', 'MarkerFaceColor', [0 0.75 0], 'MarkerEdgeColor', 'none');
    errorbar([1 2], [mean(bite_amplitudes_NI) mean(bite_amplitudes_I)],...
        [std(bite_amplitudes_NI)/sqrt(numel(bite_amplitudes_NI)) std(bite_amplitudes_I)/sqrt(numel(bite_amplitudes_I))],...
        'k', 'LineStyle', 'none');
    set(gca, 'XTick', [1 2], 'XTickLabel', {'Control', 'Inhibition'}, 'XLim', [0.4 2.6], 'TickLength', [0 0], 'FontSize', 12);
    ylabel('Normalized Sound Amplitude');
    fprintf('\n');
    try
        p = ranksum(bite_amplitudes_NI, bite_amplitudes_I, 'tail', 'both');
        disp(['bite amplitude (both): p = ' num2str(p)]);
        p = ranksum(bite_amplitudes_NI, bite_amplitudes_I, 'tail', 'right');
        disp(['bite amplitude (right): p = ' num2str(p)]);
    end
    
    if app.BiteDeviceButton.Value
        % plot bite_duration
        figure;
        hold on;
        bar(1, mean(bite_duration_NI), 0.8, 'FaceColor', [0.5 0.5 0.5]);
        plot(1+rand(1, numel(bite_duration_NI))*0.8-0.4, bite_duration_NI, 'o', 'MarkerFaceColor', [0 0 0], 'MarkerEdgeColor', 'none');
        bar(2, mean(bite_duration_I), 0.8, 'FaceColor', [0 1 0]);
        plot(2+rand(1, numel(bite_duration_I))*0.8-0.4, bite_duration_I, 'o', 'MarkerFaceColor', [0 0.75 0], 'MarkerEdgeColor', 'none');
        errorbar([1 2], [mean(bite_duration_NI) mean(bite_duration_I)],...
            [std(bite_duration_NI)/sqrt(numel(bite_duration_NI)) std(bite_duration_I)/sqrt(numel(bite_duration_I))],...
            'k', 'LineStyle', 'none');
        set(gca, 'XTick', [1 2], 'XTickLabel', {'Control', 'Inhibition'}, 'XLim', [0.4 2.6], 'TickLength', [0 0], 'FontSize', 12);
        ylabel('Total duration of all bites (s)');
        fprintf('\n');
        p = ranksum(bite_duration_NI, bite_duration_I, 'tail', 'both');
        disp(['Total duration of all bites (both): p = ' num2str(p)]);
        p = ranksum(bite_duration_NI, bite_duration_I, 'tail', 'left');
        disp(['Total duration of all bites (left): p = ' num2str(p)]);
    end
elseif app.sOnTwiceButton.Value
    bite_rate_NI = zeros(numel(trials_NI), numel(bite_edges)-1);
    bite_rate_I = zeros(numel(trials_I), numel(bite_edges)-1);
    % raster plot
    figure;
    hold on;
    for i = 1:numel(trials_NI)
        bite_rate_NI(i, :) = histcounts(trials_NI(i).bite_timestamps, bite_edges)/bite_bin;
        for j = 1:numel(trials_NI(i).bite_timestamps)
            line([trials_NI(i).bite_timestamps(j); trials_NI(i).bite_timestamps(j)],...
                [i-0.5-trials_NI(i).bite_amplitudes(j)/2; i-0.5+trials_NI(i).bite_amplitudes(j)/2],...
                'Color', 'k', 'LineStyle', '-', 'LineWidth', 1);
        end
    end
    line([4 8; 4 8], [0 0; numel(trials_NI) numel(trials_NI)], 'Color', 'r', 'LineStyle', '--', 'LineWidth', 1);
    line([13 17; 13 17], [0 0; numel(trials_NI) numel(trials_NI)], 'Color', 'r', 'LineStyle', '--', 'LineWidth', 1);
    for i = 1:numel(trials_I)
        bite_rate_I(i, :) = histcounts(trials_I(i).bite_timestamps, bite_edges)/bite_bin;
        patch([trials_I(i).laser_timestamps(1, 1) trials_I(i).laser_timestamps(1, 1) trials_I(i).laser_timestamps(2, 1) trials_I(i).laser_timestamps(2, 1)],...
            [numel(trials_NI)+i-1 numel(trials_NI)+i numel(trials_NI)+i numel(trials_NI)+i-1], [0 1 0], 'FaceAlpha', 1, 'EdgeColor', 'none');
        patch([trials_I(i).laser_timestamps(1, 2) trials_I(i).laser_timestamps(1, 2) trials_I(i).laser_timestamps(2, 2) trials_I(i).laser_timestamps(2, 2)],...
            [numel(trials_NI)+i-1 numel(trials_NI)+i numel(trials_NI)+i numel(trials_NI)+i-1], [0 1 0], 'FaceAlpha', 1, 'EdgeColor', 'none');
        for j = 1:numel(trials_I(i).bite_timestamps)
            line([trials_I(i).bite_timestamps(j); trials_I(i).bite_timestamps(j)],...
                [numel(trials_NI)+i-0.5-trials_I(i).bite_amplitudes(j)/2; numel(trials_NI)+i-0.5+trials_I(i).bite_amplitudes(j)/2],...
                'Color', 'k', 'LineStyle', '-', 'LineWidth', 1);
        end
    end
    if ~isinf(time_range(2))
        xlim(time_range);
    elseif isinf(time_range(2))
        xl = xlim;
        xl(1) = time_range(1);
        xlim(xl);
    end
    set(gca, 'YTick', 0.5:1:numel(trials_NI)+numel(trials_I)-0.5, 'YTickLabel', num2str((1:1:numel(trials_NI)+numel(trials_I))'), 'yLim', [0 numel(trials_I)+numel(trials_NI)],...
        'TickLength', [0 0], 'FontSize', 12, 'YDir', 'reverse');
    xlabel('Time (s)');
    ylabel('Trials');
    
    figure;
    plot_tj_multitrial(time_range, bite_edges(2:end)', bite_rate_NI', bite_edges(2:end)', bite_rate_I',...
        get_laser_timestamp_sound(trials_I, 'start', 1:2, 'mean'), get_laser_timestamp_sound(trials_I, 'stop', 1:2, 'mean'), 'Bite (bite/s)', '');
    
    % bar plots
    figure;
    subplot(1, 4, 1);
    hold on;
    bar(1, mean(firstafter2lastbefore_NI1), 0.8, 'FaceColor', [0.5 0.5 0.5]);
    plot(1+rand(1, numel(firstafter2lastbefore_NI1))*0.8-0.4, firstafter2lastbefore_NI1, 'o',...
        'MarkerFaceColor', [0 0 0], 'MarkerEdgeColor', 'none');
    bar(2, mean(firstafter2lastbefore_I1), 0.8, 'FaceColor', [0 1 0]);
    plot(2+rand(1, numel(firstafter2lastbefore_I1))*0.8-0.4, firstafter2lastbefore_I1, 'o',...
        'MarkerFaceColor', [0 0.75 0], 'MarkerEdgeColor', 'none');
    errorbar([1 2], [mean(firstafter2lastbefore_NI1) mean(firstafter2lastbefore_I1)],...
        [std(firstafter2lastbefore_NI1)/sqrt(numel(firstafter2lastbefore_NI1)) std(firstafter2lastbefore_I1)/sqrt(numel(firstafter2lastbefore_I1))],...
        'k', 'LineStyle', 'none');
    
    bar(3, mean(firstafter2lastbefore_NI2), 0.8, 'FaceColor', [0.5 0.5 0.5]);
    plot(3+rand(1, numel(firstafter2lastbefore_NI2))*0.8-0.4, firstafter2lastbefore_NI2, 'o',...
        'MarkerFaceColor', [0 0 0], 'MarkerEdgeColor', 'none');
    bar(4, mean(firstafter2lastbefore_I2), 0.8, 'FaceColor', [0 1 0]);
    plot(4+rand(1, numel(firstafter2lastbefore_I2))*0.8-0.4, firstafter2lastbefore_I2, 'o',...
        'MarkerFaceColor', [0 0.75 0], 'MarkerEdgeColor', 'none');
    errorbar([3 4], [mean(firstafter2lastbefore_NI2) mean(firstafter2lastbefore_I2)],...
        [std(firstafter2lastbefore_NI2)/sqrt(numel(firstafter2lastbefore_NI2)) std(firstafter2lastbefore_I2)/sqrt(numel(firstafter2lastbefore_I2))],...
        'k', 'LineStyle', 'none');
    
    set(gca, 'XTick', [1 2 3 4], 'XTickLabel', {'Control 1st', 'Inhibition 1st', 'Control 2nd', 'Inhibition 2nd'}, 'XLim', [0.4 4.6], 'TickLength', [0 0], 'FontSize', 12);
    ylabel('Interval (s)');
    fprintf('\n');
    p = ranksum(firstafter2lastbefore_NI1, firstafter2lastbefore_I1, 'tail', 'both');
    disp(['firstafter2lastbefore 1st (both): p = ' num2str(p)]);
    p = ranksum(firstafter2lastbefore_NI1, firstafter2lastbefore_I1, 'tail', 'left');
    disp(['firstafter2lastbefore 1st (left): p = ' num2str(p)]);
    p = ranksum(firstafter2lastbefore_NI2, firstafter2lastbefore_I2, 'tail', 'both');
    disp(['firstafter2lastbefore 2nd (both): p = ' num2str(p)]);
    p = ranksum(firstafter2lastbefore_NI2, firstafter2lastbefore_I2, 'tail', 'left');
    disp(['firstafter2lastbefore 2nd (left): p = ' num2str(p)]);
    
    subplot(1, 4, 2);
    hold on;
    bar(1, mean(totalbites_NI1), 0.8, 'FaceColor', [0.5 0.5 0.5]);
    plot(1+rand(1, numel(totalbites_NI1))*0.8-0.4, totalbites_NI1, 'o', 'MarkerFaceColor', [0 0 0], 'MarkerEdgeColor', 'none');
    bar(2, mean(totalbites_I1), 0.8, 'FaceColor', [0 1 0]);
    plot(2+rand(1, numel(totalbites_I1))*0.8-0.4, totalbites_I1, 'o', 'MarkerFaceColor', [0 0.75 0], 'MarkerEdgeColor', 'none');
    errorbar([1 2], [mean(totalbites_NI1) mean(totalbites_I1)],...
        [std(totalbites_NI1)/sqrt(numel(totalbites_NI1)) std(totalbites_I1)/sqrt(numel(totalbites_I1))], 'k', 'LineStyle', 'none');
    
    bar(3, mean(totalbites_NI2), 0.8, 'FaceColor', [0.5 0.5 0.5]);
    plot(3+rand(1, numel(totalbites_NI2))*0.8-0.4, totalbites_NI2, 'o', 'MarkerFaceColor', [0 0 0], 'MarkerEdgeColor', 'none');
    bar(4, mean(totalbites_I2), 0.8, 'FaceColor', [0 1 0]);
    plot(4+rand(1, numel(totalbites_I2))*0.8-0.4, totalbites_I2, 'o', 'MarkerFaceColor', [0 0.75 0], 'MarkerEdgeColor', 'none');
    errorbar([3 4], [mean(totalbites_NI2) mean(totalbites_I2)],...
        [std(totalbites_NI2)/sqrt(numel(totalbites_NI2)) std(totalbites_I2)/sqrt(numel(totalbites_I2))], 'k', 'LineStyle', 'none');
    
    set(gca, 'XTick', [1 2 3 4], 'XTickLabel', {'Control 1st', 'Inhibition 1st', 'Control 2nd', 'Inhibition 2nd'}, 'XLim', [0.4 4.6], 'TickLength', [0 0], 'FontSize', 12);
    ylabel('Bites');
    fprintf('\n');
    p = ranksum(totalbites_NI1, totalbites_I1, 'tail', 'both');
    disp(['totalbites 1st (both): p = ' num2str(p)]);
    p = ranksum(totalbites_NI1, totalbites_I1, 'tail', 'right');
    disp(['totalbites 1st (right): p = ' num2str(p)]);
    p = ranksum(totalbites_NI2, totalbites_I2, 'tail', 'both');
    disp(['totalbites 2nd (both): p = ' num2str(p)]);
    p = ranksum(totalbites_NI2, totalbites_I2, 'tail', 'right');
    disp(['totalbites 2nd (right): p = ' num2str(p)]);
    
    try
        subplot(1, 4, 3);
        hold on;
        bar(1, mean(bite_intervals_NI1), 0.8, 'FaceColor', [0.5 0.5 0.5]);
        plot(1+rand(1, numel(bite_intervals_NI1))*0.8-0.4, bite_intervals_NI1, 'o', 'MarkerFaceColor', [0 0 0], 'MarkerEdgeColor', 'none');
        bar(2, mean(bite_intervals_I1), 0.8, 'FaceColor', [0 1 0]);
        plot(2+rand(1, numel(bite_intervals_I1))*0.8-0.4, bite_intervals_I1, 'o', 'MarkerFaceColor', [0 0.75 0], 'MarkerEdgeColor', 'none');
        errorbar([1 2], [mean(bite_intervals_NI1) mean(bite_intervals_I1)],...
            [std(bite_intervals_NI1)/sqrt(numel(bite_intervals_NI1)) std(bite_intervals_I1)/sqrt(numel(bite_intervals_I1))], 'k', 'LineStyle', 'none');
        
        bar(3, mean(bite_intervals_NI2), 0.8, 'FaceColor', [0.5 0.5 0.5]);
        plot(3+rand(1, numel(bite_intervals_NI2))*0.8-0.4, bite_intervals_NI2, 'o', 'MarkerFaceColor', [0 0 0], 'MarkerEdgeColor', 'none');
        bar(4, mean(bite_intervals_I2), 0.8, 'FaceColor', [0 1 0]);
        plot(4+rand(1, numel(bite_intervals_I2))*0.8-0.4, bite_intervals_I2, 'o', 'MarkerFaceColor', [0 0.75 0], 'MarkerEdgeColor', 'none');
        errorbar([3 4], [mean(bite_intervals_NI2) mean(bite_intervals_I2)],...
            [std(bite_intervals_NI2)/sqrt(numel(bite_intervals_NI2)) std(bite_intervals_I2)/sqrt(numel(bite_intervals_I2))], 'k', 'LineStyle', 'none');
        
        set(gca, 'XTick', [1 2 3 4], 'XTickLabel', {'Control 1st', 'Inhibition 1st', 'Control 2nd', 'Inhibition 2nd'}, 'XLim', [0.4 4.6], 'TickLength', [0 0], 'FontSize', 12);
        ylabel('Bite Interval (s)');
        fprintf('\n');
        p = ranksum(bite_intervals_NI1, bite_intervals_I1, 'tail', 'both');
        disp(['bite interval 1st (both): p = ' num2str(p)]);
        p = ranksum(bite_intervals_NI1, bite_intervals_I1, 'tail', 'left');
        disp(['bite interval 1st (left): p = ' num2str(p)]);
        p = ranksum(bite_intervals_NI2, bite_intervals_I2, 'tail', 'both');
        disp(['bite interval 2nd (both): p = ' num2str(p)]);
        p = ranksum(bite_intervals_NI2, bite_intervals_I2, 'tail', 'left');
        disp(['bite interval 2nd (left): p = ' num2str(p)]);
    end
    
    try
        subplot(1, 4, 4);
        hold on;
        bar(1, mean(bite_amplitudes_NI1), 0.8, 'FaceColor', [0.5 0.5 0.5]);
        plot(1+rand(1, numel(bite_amplitudes_NI1))*0.8-0.4, bite_amplitudes_NI1, 'o', 'MarkerFaceColor', [0 0 0], 'MarkerEdgeColor', 'none');
        bar(2, mean(bite_amplitudes_I1), 0.8, 'FaceColor', [0 1 0]);
        plot(2+rand(1, numel(bite_amplitudes_I1))*0.8-0.4, bite_amplitudes_I1, 'o', 'MarkerFaceColor', [0 0.75 0], 'MarkerEdgeColor', 'none');
        errorbar([1 2], [mean(bite_amplitudes_NI1) mean(bite_amplitudes_I1)],...
            [std(bite_amplitudes_NI1)/sqrt(numel(bite_amplitudes_NI1)) std(bite_amplitudes_I1)/sqrt(numel(bite_amplitudes_I1))],...
            'k', 'LineStyle', 'none');
        
        bar(3, mean(bite_amplitudes_NI2), 0.8, 'FaceColor', [0.5 0.5 0.5]);
        plot(3+rand(1, numel(bite_amplitudes_NI2))*0.8-0.4, bite_amplitudes_NI2, 'o', 'MarkerFaceColor', [0 0 0], 'MarkerEdgeColor', 'none');
        bar(4, mean(bite_amplitudes_I2), 0.8, 'FaceColor', [0 1 0]);
        plot(4+rand(1, numel(bite_amplitudes_I2))*0.8-0.4, bite_amplitudes_I2, 'o', 'MarkerFaceColor', [0 0.75 0], 'MarkerEdgeColor', 'none');
        errorbar([3 4], [mean(bite_amplitudes_NI2) mean(bite_amplitudes_I2)],...
            [std(bite_amplitudes_NI2)/sqrt(numel(bite_amplitudes_NI2)) std(bite_amplitudes_I2)/sqrt(numel(bite_amplitudes_I2))],...
            'k', 'LineStyle', 'none');
        
        set(gca, 'XTick', [1 2 3 4], 'XTickLabel', {'Control 1st', 'Inhibition 1st', 'Control 2nd', 'Inhibition 2nd'}, 'XLim', [0.4 4.6], 'TickLength', [0 0], 'FontSize', 12);
        ylabel('Normalized Sound Amplitude');
        fprintf('\n');
        p = ranksum(bite_amplitudes_NI1, bite_amplitudes_I1, 'tail', 'both');
        disp(['bite amplitude 1st (both): p = ' num2str(p)]);
        p = ranksum(bite_amplitudes_NI1, bite_amplitudes_I1, 'tail', 'right');
        disp(['bite amplitude 1st (right): p = ' num2str(p)]);
        p = ranksum(bite_amplitudes_NI2, bite_amplitudes_I2, 'tail', 'both');
        disp(['bite amplitude 2nd (both): p = ' num2str(p)]);
        p = ranksum(bite_amplitudes_NI2, bite_amplitudes_I2, 'tail', 'right');
        disp(['bite amplitude 2nd (right): p = ' num2str(p)]);
    end
end

if app.sDelay4sOnButton.Value || app.sDelay4sOnButton_2.Value || app.WholeTrialOnButton.Value || app.BiteDeviceButton.Value
    % histograms & cdf
    try
        [~, edges] = histcounts([bite_intervals_I bite_intervals_NI]);
        figure;
        subplot(2, 2, 1);
        histogram(bite_intervals_NI, edges, 'FaceColor', [0 0 0], 'FaceAlpha', 1, 'Normalization', 'probability');
        hold on;
        histogram(bite_intervals_I, edges, 'FaceColor', 'none', 'EdgeColor', [0 1 0], 'LineWidth', 1, 'Normalization', 'probability');
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
        set(gca, 'GridLineStyle', ':', 'TickLength', [0 0], 'FontSize', 12);
        xlabel('Bite Interval (s)');
        ylabel('Cumulative Probability');
        fprintf('\n');
        if ~isempty(bite_intervals_I)
            [~, p] = kstest2(bite_intervals_NI, bite_intervals_I);
            disp(['bite_intervals (KStest): p = ' num2str(p)]);
        end
    end
    
    %                 [~, edges] = histcounts([bite_amplitudes_I bite_amplitudes_NI]);
    try
        edges = 0:0.1:1;
        subplot(2, 2, 2);
        histogram(bite_amplitudes_NI, edges, 'FaceColor', [0 0 0], 'FaceAlpha', 1, 'Normalization', 'probability');
        hold on;
        histogram(bite_amplitudes_I, edges, 'FaceColor', 'none', 'EdgeColor', [0 1 0], 'LineWidth', 1, 'Normalization', 'probability');
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
        set(gca, 'GridLineStyle', ':', 'TickLength', [0 0], 'FontSize', 12);
        xlabel('Normalized Sound Amplitude');
        ylabel('Cumulative Probability');
        fprintf('\n');
        if ~isempty(bite_amplitudes_I)
            [~, p] = kstest2(bite_amplitudes_NI, bite_amplitudes_I);
            disp(['bite_amplitudes (KStest): p = ' num2str(p)]);
        end
    end
elseif app.NoLightButton.Value
    % histograms & cdf
    [~, edges] = histcounts(bite_intervals);
    figure;
    subplot(2, 2, 1);
    histogram(bite_intervals, edges, 'FaceColor', [0 0 0], 'FaceAlpha', 1, 'Normalization', 'probability');
    set(gca, 'TickLength', [0 0], 'FontSize', 12);
    box off;
    xlabel('Bite Interval (s)');
    ylabel('Probability');
    subplot(2, 2, 3);
    hold on;
    h = cdfplot(bite_intervals);
    set(h, 'Color', [0 0 0], 'LineWidth', 1.5);
    set(gca, 'GridLineStyle', ':', 'TickLength', [0 0], 'FontSize', 12);
    xlabel('Bite Interval (s)');
    ylabel('Cumulative Probability');
    
    %                 [~, edges] = histcounts(bite_amplitudes);
    edges = 0:0.1:1;
    subplot(2, 2, 2);
    histogram(bite_amplitudes, edges, 'FaceColor', [0 0 0], 'FaceAlpha', 1, 'Normalization', 'probability');
    set(gca, 'TickLength', [0 0], 'FontSize', 12);
    box off;
    xlabel('Normalized Sound Amplitude');
    ylabel('Probability');
    subplot(2, 2, 4);
    hold on;
    h = cdfplot(bite_amplitudes);
    set(h, 'Color', [0 0 0], 'LineWidth', 1.5);
    set(gca, 'GridLineStyle', ':', 'TickLength', [0 0], 'FontSize', 12);
    xlabel('Normalized Sound Amplitude');
    ylabel('Cumulative Probability');
elseif app.sOnTwiceButton.Value
    % histograms & cdf
    try
        [~, edges] = histcounts([bite_intervals_I1 bite_intervals_NI1 bite_intervals_I2 bite_intervals_NI2]);
        figure;
        subplot(3, 2, 1);
        histogram(bite_intervals_NI1, edges, 'FaceColor', [0 0 0], 'FaceAlpha', 1, 'Normalization', 'probability');
        hold on;
        histogram(bite_intervals_I1, edges, 'FaceColor', 'none', 'EdgeColor', [0 1 0], 'LineWidth', 1, 'Normalization', 'probability');
        set(gca, 'TickLength', [0 0], 'FontSize', 12);
        box off;
        xlabel('Bite Interval (s)');
        ylabel('Probability');
        title('First 4 seconds');
        
        subplot(3, 2, 3);
        histogram(bite_intervals_NI2, edges, 'FaceColor', [0 0 0], 'FaceAlpha', 1, 'Normalization', 'probability');
        hold on;
        histogram(bite_intervals_I2, edges, 'FaceColor', 'none', 'EdgeColor', [0 1 0], 'LineWidth', 1, 'Normalization', 'probability');
        set(gca, 'TickLength', [0 0], 'FontSize', 12);
        box off;
        xlabel('Bite Interval (s)');
        ylabel('Probability');
        title('Second 4 seconds');
        
        subplot(3, 2, 5);
        hold on;
        h1 = cdfplot(bite_intervals_NI1);
        set(h1, 'Color', [0.75 0.75 0.75], 'LineWidth', 1.5);
        if ~isempty(bite_intervals_I1)
            h2 = cdfplot(bite_intervals_I1);
            set(h2, 'Color', [0 1 0], 'LineWidth', 1.5);
        end
        h3 = cdfplot(bite_intervals_NI2);
        set(h3, 'Color', [0 0 0], 'LineWidth', 1.5);
        if ~isempty(bite_intervals_I2)
            h4 = cdfplot(bite_intervals_I2);
            set(h4, 'Color', [0 0.5 0], 'LineWidth', 1.5);
        end
        set(gca, 'GridLineStyle', ':', 'TickLength', [0 0], 'FontSize', 12);
        xlabel('Bite Interval (s)');
        ylabel('Cumulative Probability');
        legend([h1 h2 h3 h4], {'Control 1st', 'Inhibition 1st', 'Control 2nd', 'Inhibition 2nd'}, 'FontSize', 12);
        legend('boxoff');
        fprintf('\n');
        if ~isempty(bite_intervals_I1)
            [~, p] = kstest2(bite_intervals_NI1, bite_intervals_I1);
            disp(['bite_intervals 1st (KStest): p = ' num2str(p)]);
        end
        if ~isempty(bite_intervals_I2)
            [~, p] = kstest2(bite_intervals_NI2, bite_intervals_I2);
            disp(['bite_intervals 2nd (KStest): p = ' num2str(p)]);
        end
    end
    
    %                 [~, edges] = histcounts([bite_amplitudes_I bite_amplitudes_NI]);
    try
        edges = 0:0.1:1;
        subplot(3, 2, 2);
        histogram(bite_amplitudes_NI1, edges, 'FaceColor', [0 0 0], 'FaceAlpha', 1, 'Normalization', 'probability');
        hold on;
        histogram(bite_amplitudes_I1, edges, 'FaceColor', 'none', 'EdgeColor', [0 1 0], 'LineWidth', 1, 'Normalization', 'probability');
        set(gca, 'TickLength', [0 0], 'FontSize', 12);
        box off;
        xlabel('Normalized Sound Amplitude');
        ylabel('Probability');
        title('First 4 seconds');
        
        subplot(3, 2, 4);
        histogram(bite_amplitudes_NI2, edges, 'FaceColor', [0 0 0], 'FaceAlpha', 1, 'Normalization', 'probability');
        hold on;
        histogram(bite_amplitudes_I2, edges, 'FaceColor', 'none', 'EdgeColor', [0 1 0], 'LineWidth', 1, 'Normalization', 'probability');
        set(gca, 'TickLength', [0 0], 'FontSize', 12);
        box off;
        xlabel('Normalized Sound Amplitude');
        ylabel('Probability');
        title('Second 4 seconds');
        
        subplot(3, 2, 6);
        hold on;
        h1 = cdfplot(bite_amplitudes_NI1);
        set(h1, 'Color', [0.75 0.75 0.75], 'LineWidth', 1.5);
        if ~isempty(bite_amplitudes_I1)
            h2 = cdfplot(bite_amplitudes_I1);
            set(h2, 'Color', [0 1 0], 'LineWidth', 1.5);
        end
        h3 = cdfplot(bite_amplitudes_NI2);
        set(h3, 'Color', [0 0 0], 'LineWidth', 1.5);
        if ~isempty(bite_amplitudes_I2)
            h4 = cdfplot(bite_amplitudes_I2);
            set(h4, 'Color', [0 0.5 0], 'LineWidth', 1.5);
        end
        set(gca, 'GridLineStyle', ':', 'TickLength', [0 0], 'FontSize', 12);
        xlabel('Normalized Sound Amplitude');
        ylabel('Cumulative Probability');
        legend([h1 h2 h3 h4], {'Control 1st', 'Inhibition 1st', 'Control 2nd', 'Inhibition 2nd'}, 'FontSize', 12);
        legend('boxoff');
        fprintf('\n');
        if ~isempty(bite_amplitudes_I1)
            [~, p] = kstest2(bite_amplitudes_NI1, bite_amplitudes_I1);
            disp(['bite_amplitudes 1st (KStest): p = ' num2str(p)]);
        end
        if ~isempty(bite_amplitudes_I2)
            [~, p] = kstest2(bite_amplitudes_NI2, bite_amplitudes_I2);
            disp(['bite_amplitudes 2nd (KStest): p = ' num2str(p)]);
        end
    end
end

if app.sDelay4sOnButton.Value
    result.firstafter2lastbefore_I = firstafter2lastbefore_I;
    result.firstafter2lastbefore_NI = firstafter2lastbefore_NI;
    result.totalbites_I = totalbites_I;
    result.totalbites_NI = totalbites_NI;
    result.bite_intervals_I = bite_intervals_I;
    result.bite_intervals_NI = bite_intervals_NI;
    result.bite_amplitudes_I =bite_amplitudes_I;
    result.bite_amplitudes_NI =bite_amplitudes_NI;
    result.trials_I = trials_I;
    result.trials_NI = trials_NI;
elseif app.sDelay4sOnButton_2.Value || app.WholeTrialOnButton.Value || app.BiteDeviceButton.Value
    result.firstbite_time_I = firstbite_time_I;
    result.firstbite_time_NI = firstbite_time_NI;
    result.totalbites_I = totalbites_I;
    result.totalbites_NI = totalbites_NI;
    if app.BiteDeviceButton.Value
        result.bite_duration_I = bite_duration_I;
        result.bite_duration_NI = bite_duration_NI;
    end
    result.bite_intervals_I = bite_intervals_I;
    result.bite_intervals_NI = bite_intervals_NI;
    result.bite_amplitudes_I =bite_amplitudes_I;
    result.bite_amplitudes_NI =bite_amplitudes_NI;
    result.trials_I = trials_I;
    result.trials_NI = trials_NI;
elseif app.NoLightButton.Value
    result.firstbite_time = firstbite_time;
    result.all_trials = all_trials;
    result.totalbites = totalbites;
    result.bite_intervals = bite_intervals;
    result.bite_amplitudes = bite_amplitudes;
elseif app.sOnTwiceButton.Value
    result.firstafter2lastbefore_I1 = firstafter2lastbefore_I1;
    result.firstafter2lastbefore_NI1 = firstafter2lastbefore_NI1;
    result.totalbites_I1 = totalbites_I1;
    result.totalbites_NI1 = totalbites_NI1;
    result.bite_intervals_I1 = bite_intervals_I1;
    result.bite_intervals_NI1 = bite_intervals_NI1;
    result.bite_amplitudes_I1 =bite_amplitudes_I1;
    result.bite_amplitudes_NI1 =bite_amplitudes_NI1;
    result.firstafter2lastbefore_I2 = firstafter2lastbefore_I2;
    result.firstafter2lastbefore_NI2 = firstafter2lastbefore_NI2;
    result.totalbites_I2 = totalbites_I2;
    result.totalbites_NI2 = totalbites_NI2;
    result.bite_intervals_I2 = bite_intervals_I2;
    result.bite_intervals_NI2 = bite_intervals_NI2;
    result.bite_amplitudes_I2 =bite_amplitudes_I2;
    result.bite_amplitudes_NI2 =bite_amplitudes_NI2;
    result.trials_I = trials_I;
    result.trials_NI = trials_NI;
end
[file, path] = uiputfile('*.mat', 'Save result');
if file ~= 0
    save([path file], 'result');
end