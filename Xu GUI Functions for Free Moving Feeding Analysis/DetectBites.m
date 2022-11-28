function DetectBites(app, Exp_Path, SampleRate, lowpassfreq, highpassfreq, Video_annotation, standard_delay)
minbite_interval = 0.005;
trial = app.TrialsListBox.Value;
audiolocation = Exp_Path(1:end-7);
[blow, alow] = butter(6, lowpassfreq/(SampleRate/2), 'low');
[bhigh, ahigh] = butter(6, highpassfreq/(SampleRate/2), 'high');
try
    temp = load([audiolocation '\Detected_Bite_Events.mat']);
    Audio_analysis = temp.Audio_analysis;
catch
    Audio_analysis = struct('time_bites', cell(app.TrialsListBox.ItemsData(end), 1), 'amplitude_bites',...
        cell(app.TrialsListBox.ItemsData(end), 1));
    save([audiolocation '\Detected_Bite_Events.mat'], 'Audio_analysis');
end
workbar(0, 'Computing Ongoing...', 'Progress');
if app.CheckDetecetionCheckBox.Value
    for i = 1:numel(trial)
        if ~Video_annotation(trial(i)).Disgard
            try
                if app.BiteDeviceButton.Value
                    temp = load([Exp_Path '\LabelledEvents' num2str(trial(i)) '.mat']);
                    LabelledEvents = temp.LabelledEvents;
                end
                
                temp = load([audiolocation '\Trial' num2str(trial(i), '%03d')]);
                data = temp.audiodata;
                figure;
                subplot(3, 1, 1);
                if isfield(data, 'laser_starttime') && app.OptSessionCheckBox.Value
                    if ~isempty(data.laser_starttime)
                        for j = 1:numel(data.laser_starttime)
                            patch([data.laser_starttime(j) data.laser_starttime(j) data.laser_endtime(j) data.laser_endtime(j)],...
                                [min(data.signal) max(data.signal) max(data.signal) min(data.signal)],...
                                [0 1 0], 'FaceAlpha', 1, 'EdgeColor', 'none');
                        end
                    else
                        line([standard_delay; standard_delay], [min(data.signal); max(data.signal)],...
                            'Color', 'g', 'LineStyle', '--', 'LineWidth', 1);
                    end
                end
                hold on;
                plot(data.timestamps, data.signal, '-k');
                set(gca, 'xLim', [0 data.timestamps(end)], 'yLim', [min(data.signal) max(data.signal)], 'TickDir', 'out');
                box off;
                xlabel('Time (s)');
                ylabel('Voltage (V)');
                title(['Trial ' num2str(trial(i)) ' (raw signal)']);
                
                data.signal = filtfilt(blow, alow, data.signal);
                data.signal = filtfilt(bhigh, ahigh, data.signal);
                subplot(3, 1, 2);
                if isfield(data, 'laser_starttime') && app.OptSessionCheckBox.Value
                    if ~isempty(data.laser_starttime)
                        for j = 1:numel(data.laser_starttime)
                            patch([data.laser_starttime(j) data.laser_starttime(j) data.laser_endtime(j) data.laser_endtime(j)],...
                                [min(data.signal) max(data.signal) max(data.signal) min(data.signal)],...
                                [0 1 0], 'FaceAlpha', 1, 'EdgeColor', 'none');
                        end
                    else
                        line([standard_delay; standard_delay], [min(data.signal); max(data.signal)],...
                            'Color', 'g', 'LineStyle', '--', 'LineWidth', 1);
                    end
                end
                hold on;
                plot(data.timestamps, data.signal, '-k');
                temp = abs(data.signal);
                plot(data.timestamps, smoothdata(temp, 'gaussian', minbite_interval*SampleRate), '-r');
                set(gca, 'xLim', [0 data.timestamps(end)], 'yLim', [min(data.signal) max(data.signal)], 'TickDir', 'out');
                box off;
                xlabel('Time (s)');
                ylabel('Voltage (V)');
                title(['Trial ' num2str(trial(i)) ' (filtered signal)']);
                
                data.signal = smoothdata(temp, 'gaussian', minbite_interval*SampleRate);
                subplot(3, 1, 3);
                if isfield(data, 'laser_starttime') && app.OptSessionCheckBox.Value
                    if ~isempty(data.laser_starttime)
                        for j = 1:numel(data.laser_starttime)
                            patch([data.laser_starttime(j) data.laser_starttime(j) data.laser_endtime(j) data.laser_endtime(j)],...
                                [min(data.signal) max(data.signal) max(data.signal) min(data.signal)],...
                                [0 1 0], 'FaceAlpha', 1, 'EdgeColor', 'none');
                        end
                    else
                        line([standard_delay; standard_delay], [min(data.signal); max(data.signal)],...
                            'Color', 'g', 'LineStyle', '--', 'LineWidth', 1);
                    end
                end
                hold on;
                plot(data.timestamps, data.signal, '-k');
                threshold = mean(data.signal)+app.ThresholdEditField.Value*std(data.signal);
                line([0; data.timestamps(end)], [threshold; threshold], 'Color', 'r', 'LineStyle', '-');
                if ~app.BiteDeviceButton.Value
                    if ~isempty(Video_annotation(trial(i)).time_ready2bite)
                        line([Video_annotation(trial(i)).time_ready2bite; Video_annotation(trial(i)).time_ready2bite],...
                            [min(data.signal); max(data.signal)], 'Color', 'k', 'LineStyle', '--');
                    end
                    if ~isempty(Video_annotation(trial(i)).time_feeding_end)
                        line([Video_annotation(trial(i)).time_feeding_end; Video_annotation(trial(i)).time_feeding_end],...
                            [min(data.signal); max(data.signal)], 'Color', 'k', 'LineStyle', '--');
                    end
                else
                    line([LabelledEvents.MouthRetrievalStart'; LabelledEvents.MouthRetrievalStart'],...
                        repmat([min(data.signal); max(data.signal)], 1, numel(LabelledEvents.MouthRetrievalStart)), 'Color', 'k', 'LineStyle', '--');
                    line([LabelledEvents.MouthRetrievalEnd'; LabelledEvents.MouthRetrievalEnd'],...
                        repmat([min(data.signal); max(data.signal)], 1, numel(LabelledEvents.MouthRetrievalEnd)), 'Color', 'k', 'LineStyle', '--');
                end
                events = find(data.signal >= threshold);
                bite_ids = [];
                n = 1;
                bite_amplitude = [];
                bite_time = [];
                for j = 1:numel(events)
                    if isempty(bite_ids)
                        bite_ids = events(j);
                    else
                        if events(j)-bite_ids(end) == 1
                            bite_ids = [bite_ids events(j)];
                            if j == numel(events)
                                bite_signal = data.signal(bite_ids);
                                bite_timestamps = data.timestamps(bite_ids);
                                [bite_amplitude(n), index] = max(bite_signal);
                                bite_time(n) = bite_timestamps(index);
                            end
                        else
                            bite_signal = data.signal(bite_ids);
                            bite_timestamps = data.timestamps(bite_ids);
                            [bite_amplitude(n), index] = max(bite_signal);
                            bite_time(n) = bite_timestamps(index);
                            n = n+1;
                            bite_ids = events(j);
                        end
                    end
                end
                if ~app.BiteDeviceButton.Value
                    if ~isempty(Video_annotation(trial(i)).time_ready2bite) && ~isempty(Video_annotation(trial(i)).time_feeding_end)
                        bite_amplitude(bite_time < Video_annotation(trial(i)).time_ready2bite |...
                            bite_time > Video_annotation(trial(i)).time_feeding_end) = [];
                        bite_time(bite_time < Video_annotation(trial(i)).time_ready2bite |...
                            bite_time > Video_annotation(trial(i)).time_feeding_end) = [];
                    elseif ~isempty(Video_annotation(trial(i)).time_ready2bite)
                        bite_amplitude(bite_time < Video_annotation(trial(i)).time_ready2bite) = [];
                        bite_time(bite_time < Video_annotation(trial(i)).time_ready2bite) = [];
                    elseif ~isempty(Video_annotation(trial(i)).time_feeding_end)
                        bite_amplitude(bite_time > Video_annotation(trial(i)).time_feeding_end) = [];
                        bite_time(bite_time > Video_annotation(trial(i)).time_feeding_end) = [];
                    end
                else
                    bite_amplitude_temp = [];
                    bite_time_temp = [];
                    for j = 1:numel(LabelledEvents.MouthRetrievalStart)
                        bite_amplitude_temp = [bite_amplitude_temp bite_amplitude(bite_time < LabelledEvents.MouthRetrievalEnd(j) & bite_time > LabelledEvents.MouthRetrievalStart(j))];
                        bite_time_temp = [bite_time_temp bite_time(bite_time < LabelledEvents.MouthRetrievalEnd(j) & bite_time > LabelledEvents.MouthRetrievalStart(j))];
                    end
                    bite_amplitude = bite_amplitude_temp;
                    bite_time = bite_time_temp;
                end
                plot(bite_time, bite_amplitude, 'ob');
                set(gca, 'xLim', [0 data.timestamps(end)], 'yLim', [min(data.signal) max(data.signal)], 'TickDir', 'out');
                box off;
                xlabel('Time (s)');
                ylabel('Voltage (V)');
                title(['Trial ' num2str(trial(i)) ' (filtered and smoothed signal)']);
                
                if isempty(bite_time)
                    Audio_analysis(trial(i)).time_bites = [];
                    Audio_analysis(trial(i)).amplitude_bites = [];
                else
                    Audio_analysis(trial(i)).time_bites = bite_time;
                    Audio_analysis(trial(i)).amplitude_bites = bite_amplitude;
                end
                if isfield(data, 'laser_starttime')
                    Audio_analysis(trial(i)).laser_timestamps = [data.laser_starttime; data.laser_endtime];
                end
                if app.BiteDeviceButton.Value
                    Audio_analysis(trial(i)).bite_window = [LabelledEvents.MouthRetrievalStart LabelledEvents.MouthRetrievalEnd];
                end
                Audio_analysis(trial(i)).threshold = app.ThresholdEditField.Value;
            end
        end
        workbar(i/numel(trial), [num2str(i) '/' num2str(numel(trial))], 'Progress');
    end
else
    for i = 1:numel(trial)
        if ~Video_annotation(trial(i)).Disgard
            try
                if app.BiteDeviceButton.Value
                    temp = load([Exp_Path '\LabelledEvents' num2str(trial(i)) '.mat']);
                    LabelledEvents = temp.LabelledEvents;
                end
                
                temp = load([audiolocation '\Trial' num2str(trial(i), '%03d')]);
                data = temp.audiodata;
                data.signal = filtfilt(blow, alow, data.signal);
                data.signal = filtfilt(bhigh, ahigh, data.signal);
                data.signal = abs(data.signal);
                data.signal = smoothdata(data.signal, 'gaussian', minbite_interval*SampleRate);
                
                threshold = mean(data.signal)+app.ThresholdEditField.Value*std(data.signal);
                events = find(data.signal >= threshold);
                bite_ids = [];
                n = 1;
                bite_amplitude = [];
                bite_time = [];
                for j = 1:numel(events)
                    if isempty(bite_ids)
                        bite_ids = events(j);
                    else
                        if events(j)-bite_ids(end) == 1
                            bite_ids = [bite_ids events(j)];
                            if j == numel(events)
                                bite_signal = data.signal(bite_ids);
                                bite_timestamps = data.timestamps(bite_ids);
                                [bite_amplitude(n), index] = max(bite_signal);
                                bite_time(n) = bite_timestamps(index);
                            end
                        else
                            bite_signal = data.signal(bite_ids);
                            bite_timestamps = data.timestamps(bite_ids);
                            [bite_amplitude(n), index] = max(bite_signal);
                            bite_time(n) = bite_timestamps(index);
                            n = n+1;
                            bite_ids = events(j);
                        end
                    end
                end
                if ~app.BiteDeviceButton.Value
                    if ~isempty(Video_annotation(trial(i)).time_ready2bite) && ~isempty(Video_annotation(trial(i)).time_feeding_end)
                        bite_amplitude(bite_time < Video_annotation(trial(i)).time_ready2bite |...
                            bite_time > Video_annotation(trial(i)).time_feeding_end) = [];
                        bite_time(bite_time < Video_annotation(trial(i)).time_ready2bite |...
                            bite_time > Video_annotation(trial(i)).time_feeding_end) = [];
                    elseif ~isempty(Video_annotation(trial(i)).time_ready2bite)
                        bite_amplitude(bite_time < Video_annotation(trial(i)).time_ready2bite) = [];
                        bite_time(bite_time < Video_annotation(trial(i)).time_ready2bite) = [];
                    elseif ~isempty(Video_annotation(trial(i)).time_feeding_end)
                        bite_amplitude(bite_time > Video_annotation(trial(i)).time_feeding_end) = [];
                        bite_time(bite_time > Video_annotation(trial(i)).time_feeding_end) = [];
                    end
                else
                    bite_amplitude_temp = [];
                    bite_time_temp = [];
                    for j = 1:numel(LabelledEvents.MouthRetrievalStart)
                        bite_amplitude_temp = [bite_amplitude_temp bite_amplitude(bite_time < LabelledEvents.MouthRetrievalEnd(j) & bite_time > LabelledEvents.MouthRetrievalStart(j))];
                        bite_time_temp = [bite_time_temp bite_time(bite_time < LabelledEvents.MouthRetrievalEnd(j) & bite_time > LabelledEvents.MouthRetrievalStart(j))];
                    end
                    bite_amplitude = bite_amplitude_temp;
                    bite_time = bite_time_temp;
                end
                if isempty(bite_time)
                    Audio_analysis(trial(i)).time_bites = [];
                    Audio_analysis(trial(i)).amplitude_bites = [];
                else
                    Audio_analysis(trial(i)).time_bites = bite_time;
                    Audio_analysis(trial(i)).amplitude_bites = bite_amplitude;
                end
                if isfield(data, 'laser_starttime')
                    Audio_analysis(trial(i)).laser_timestamps = [data.laser_starttime; data.laser_endtime];
                end
                if app.BiteDeviceButton.Value
                    Audio_analysis(trial(i)).bite_window = [LabelledEvents.MouthRetrievalStart LabelledEvents.MouthRetrievalEnd];
                end
                Audio_analysis(trial(i)).threshold = app.ThresholdEditField.Value;
            end
        end
        workbar(i/numel(trial), [num2str(i) '/' num2str(numel(trial))], 'Progress');
    end
end
save([audiolocation '\Detected_Bite_Events.mat'], 'Audio_analysis');
msgbox('Done !');