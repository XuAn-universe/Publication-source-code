function EventPlot(app, Exp_Path, FrameRate, nmedian)
value = app.TrialsListBox.Value;
value = sort(value);
trials = numel(value);
maxtime = zeros(4, trials);
trials_I = [];
trials_NI = [];
trial_category = zeros(2, trials);
laser_timestamps = [];
njaw4retrieval1 = nan(1, trials);
RetrievalStart1st = nan(1, trials);

% get bite events
Bite_events = [];
try
    audiolocation = Exp_Path(1:end-7);
    temp = load([audiolocation '\Detected_Bite_Events.mat']);
    Bite_events = temp.Audio_analysis;
end

% process trajectories
if app.TrackingDataCheckBox.Value
    label_table = table2array(app.UITable.Data);
else
    helpdlg('''Tracking Data'' needs to be checked to display trajectories');
end
dz_max = -inf;
dz_min = inf;
z_max = -inf;
z_min = inf;

% get photometry data
fpdata_all = [];
try
    fpdata_all = load([audiolocation '\FPData.mat']);
    nchannel = size(fpdata_all.zsignal_all, 1);
    for i = 1:nchannel
        lgdtext{i} = ['Channel ' num2str(i)];
    end
end
zscore_max = -inf;
zscore_min = inf;

hsp = cell(1, 4);
figure;
for i = 1:4
    hsp{i} = subplot(4, 1, i);
    hold on;
end
for i = 1:trials
    LabelledEvents = [];
    try
        temp = load([Exp_Path '\LabelledEvents' num2str(value(i)) '.mat']);
        LabelledEvents = temp.LabelledEvents;
%         LabelledEvents.FeedingEnd = LabelledEvents.FeedingEnd+1; %
        if ~app.OptSessionCheckBox.Value
            MouthRetrievalStart = LabelledEvents.MouthRetrievalStart;
            RetrievalStart = LabelledEvents.RetrievalStart;
            SitEnd = LabelledEvents.SitEnd;
            maxtime(1, i) = LabelledEvents.FeedingEnd;
        elseif app.sDelay4sOnButton.Value
            maxtime(1, i) = max(LabelledEvents.BiteBoutStart);
        elseif app.sDelay4sOnButton_2.Value
            maxtime(1, i) = 0;
        end
    end
    
    bite_timestamps = [];
    if ~isempty(Bite_events)
        bite_timestamps = Bite_events(value(i)).time_bites;
        if ~isempty(bite_timestamps)
%             bite_amplitudes = Bite_events(value(i)).amplitude_bites/(max(Bite_events(value(i)).amplitude_bites));
            bite_amplitudes = Bite_events(value(i)).amplitude_bites./Bite_events(value(i)).amplitude_bites;
            if app.OptSessionCheckBox.Value
                if app.sDelay4sOnButton.Value
                    bite_amplitudes(bite_timestamps > max(LabelledEvents.BiteBoutStart)) = [];
                    bite_timestamps(bite_timestamps > max(LabelledEvents.BiteBoutStart)) = [];
                elseif app.sDelay4sOnButton_2.Value
                    bite_amplitudes(bite_timestamps > bite_timestamps(1)) = [];
                    bite_timestamps(bite_timestamps > bite_timestamps(1)) = [];
                end
            end
            maxtime(2, i) = max(bite_timestamps);
        end
    end
    
    if app.TrackingDataCheckBox.Value
        try
            [~, ~, ~, ~, ~, x_nose, y_nose, z_nose, speed_nose, acceleration_nose, laserstart, laserstop] = trajectory_postprocessing(3, Exp_Path, value(i),...
                label_table, nmedian, FrameRate);
            [~, ~, ~, ~, ~, x_pawl, y_pawl, z_pawl, speed_pawl, acceleration_pawl, ~, ~] = trajectory_postprocessing(10, Exp_Path, value(i), label_table, nmedian, FrameRate);
            [~, ~, ~, ~, ~, x_pawr, y_pawr, z_pawr, speed_pawr, acceleration_pawr, ~, ~] = trajectory_postprocessing(16, Exp_Path, value(i), label_table, nmedian, FrameRate);
            
            pawl2nose = z_pawl(:, 1)-z_nose(:, 1);
            pawr2nose = z_pawr(:, 1)-z_nose(:, 1);
            t_trajectory = (1:size(z_nose, 1))'/FrameRate;
            if ~app.OptSessionCheckBox.Value
                if ~isempty(LabelledEvents.FeedingEnd)
                    pawl2nose(t_trajectory > LabelledEvents.FeedingEnd) = [];
                    pawr2nose(t_trajectory > LabelledEvents.FeedingEnd) = [];
                    z_nose(t_trajectory > LabelledEvents.FeedingEnd, :) = [];
                    z_pawl(t_trajectory > LabelledEvents.FeedingEnd, :) = [];
                    z_pawr(t_trajectory > LabelledEvents.FeedingEnd, :) = [];
                elseif ~isempty(bite_timestamps)
                    pawl2nose(t_trajectory > max(bite_timestamps)) = [];
                    pawr2nose(t_trajectory > max(bite_timestamps)) = [];
                    z_nose(t_trajectory > max(bite_timestamps), :) = [];
                    z_pawl(t_trajectory > max(bite_timestamps), :) = [];
                    z_pawr(t_trajectory > max(bite_timestamps), :) = [];
                end
            elseif app.sDelay4sOnButton.Value
                pawl2nose(t_trajectory > max(LabelledEvents.BiteBoutStart)) = [];
                pawr2nose(t_trajectory > max(LabelledEvents.BiteBoutStart)) = [];
                z_nose(t_trajectory > max(LabelledEvents.BiteBoutStart), :) = [];
                z_pawl(t_trajectory > max(LabelledEvents.BiteBoutStart), :) = [];
                z_pawr(t_trajectory > max(LabelledEvents.BiteBoutStart), :) = [];
            elseif app.sDelay4sOnButton_2.Value
                if ~isempty(bite_timestamps)
                    pawl2nose(t_trajectory > max(bite_timestamps)) = [];
                    pawr2nose(t_trajectory > max(bite_timestamps)) = [];
                    z_nose(t_trajectory > max(bite_timestamps), :) = [];
                    z_pawl(t_trajectory > max(bite_timestamps), :) = [];
                    z_pawr(t_trajectory > max(bite_timestamps), :) = [];
                end
            end
            dz_max = max(dz_max, max([pawl2nose; pawr2nose]));
            dz_min = min(dz_min, min([pawl2nose; pawr2nose]));
            z_max = max(z_max, max([z_nose(:, 1); z_pawl(:, 1); z_pawr(:, 1)]));
            z_min = min(z_min, min([z_nose(:, 1); z_pawl(:, 1); z_pawr(:, 1)]));
        end
    end
    
    if ~isempty(fpdata_all)
        fpdata = fpdata_all.zsignal_all(:, value(i));
        fpdata_zsignal = [];
        fpdata_t = [];
        for j = 1:nchannel
            fpdata_zsignal(:, j) = fpdata{j}(:, 2);
        end
        fpdata_t = fpdata{1}(:, 1);
        if ~isempty(LabelledEvents.FeedingEnd)
            fpdata_zsignal(fpdata_t > LabelledEvents.FeedingEnd(end)+4, :) = [];
        elseif ~isempty(bite_timestamps)
            fpdata_zsignal(fpdata_t > max(bite_timestamps)+4, :) = [];
        end
        zscore_max = max(zscore_max, max(fpdata_zsignal(:)));
        zscore_min = min(zscore_min, min(fpdata_zsignal(:)));
    end
    
    if ~app.OptSessionCheckBox.Value
        for ii = 1:4
            if ~isempty(LabelledEvents)
                plot_events(hsp{ii}, i, trials, LabelledEvents);
            end
            
            if ~isempty(bite_timestamps)
                for j = 1:numel(bite_timestamps)
                    line(hsp{ii}, [bite_timestamps(j) bite_timestamps(j)], [(trials-i)*6+0.5-bite_amplitudes(j)/2 (trials-i)*6+0.5+bite_amplitudes(j)/2],...
                        'Color', [1 0 1], 'LineStyle', '-', 'LineWidth', 1);
                end
            end
        end
        
        if ~isempty(MouthRetrievalStart) && ~isempty(RetrievalStart) && ~isempty(SitEnd)
            if any(MouthRetrievalStart == RetrievalStart(1))
                njaw4retrieval1(i) = sum(MouthRetrievalStart < SitEnd(1));
            end
        end
        if ~isempty(RetrievalStart) && ~isempty(SitEnd)
            if RetrievalStart(1) < SitEnd(1)
                RetrievalStart1st(i) = RetrievalStart(1);
            end
        end
    else
        if isempty(Bite_events(value(i)).laser_timestamps)
            trials_NI(end+1).bite_timestamps = bite_timestamps;
            if ~isempty(bite_timestamps)
                trials_NI(end).bite_amplitudes = bite_amplitudes;
            end
            trials_NI(end).LabelledEvents = LabelledEvents;
            trial_category(2, i) = numel(trials_NI);
        else
            trials_I(end+1).bite_timestamps = bite_timestamps;
            if ~isempty(bite_timestamps)
                trials_I(end).bite_amplitudes = bite_amplitudes;
            end
            trials_I(end).LabelledEvents = LabelledEvents;
            trial_category(1, i) = 1;
            trial_category(2, i) = numel(trials_I);
            laser_timestamps = Bite_events(value(i)).laser_timestamps;
        end
    end
end

try
    assignin('base', 'dz_range', dz_max-dz_min);
    assignin('base', 'z_range', z_max-z_min);
end
if ~isempty(fpdata_all)
    assignin('base', 'zscore_range', zscore_max-zscore_min);
    assignin('base', 'zscore_max', zscore_max);
    assignin('base', 'zscore_min', zscore_min);
end

if app.OptSessionCheckBox.Value
    assignin('base', 'trials_NI', trials_NI);
    assignin('base', 'trials_I', trials_I);
end

if ~app.OptSessionCheckBox.Value
    result.njaw4retrieval1 = njaw4retrieval1;
    result.RetrievalStart1st = RetrievalStart1st;
    assignin('base', 'result', result);
end

if dz_max ~= -inf && dz_min ~= inf
    z2y = 6/(z_max-z_min);
    dz2y = 6/(dz_max-dz_min);
    for i = 1:trials
        if app.TrackingDataCheckBox.Value
            try
                [~, ~, ~, ~, ~, x_nose, y_nose, z_nose, speed_nose, acceleration_nose, laserstart, laserstop] = trajectory_postprocessing(3, Exp_Path, value(i),...
                    label_table, nmedian, FrameRate);
                [~, ~, ~, ~, ~, x_pawl, y_pawl, z_pawl, speed_pawl, acceleration_pawl, ~, ~] = trajectory_postprocessing(10, Exp_Path, value(i), label_table, nmedian, FrameRate);
                [~, ~, ~, ~, ~, x_pawr, y_pawr, z_pawr, speed_pawr, acceleration_pawr, ~, ~] = trajectory_postprocessing(16, Exp_Path, value(i), label_table, nmedian, FrameRate);
%                 [~, ~, ~, z_pawr, ~, x_pawr, y_pawr, ~, speed_pawr, acceleration_pawr, ~, ~] = trajectory_postprocessing(16, Exp_Path, value(i), label_table, nmedian, FrameRate); %
                
                pawl2nose = z_pawl(:, 1)-z_nose(:, 1);
                pawr2nose = z_pawr(:, 1)-z_nose(:, 1);
                t_trajectory = (1:size(z_nose, 1))'/FrameRate;
                if maxtime(1, i) ~= 0
                    pawl2nose(t_trajectory > maxtime(1, i)) = [];
                    pawr2nose(t_trajectory > maxtime(1, i)) = [];
                    z_nose(t_trajectory > maxtime(1, i), :) = [];
                    z_pawl(t_trajectory > maxtime(1, i), :) = [];
                    z_pawr(t_trajectory > maxtime(1, i), :) = [];
                    t_trajectory(t_trajectory > maxtime(1, i)) = [];
                elseif maxtime(2, i) ~= 0
                    pawl2nose(t_trajectory > maxtime(2, i)) = [];
                    pawr2nose(t_trajectory > maxtime(2, i)) = [];
                    z_nose(t_trajectory > maxtime(2, i), :) = [];
                    z_pawl(t_trajectory > maxtime(2, i), :) = [];
                    z_pawr(t_trajectory > maxtime(2, i), :) = [];
                    t_trajectory(t_trajectory > maxtime(2, i)) = [];
                end
                maxtime(3, i) = t_trajectory(end);
                if ~app.OptSessionCheckBox.Value
                    temp = (trials-i)*6;
                else
                    if trial_category(1, i)
                        temp = (trials-numel(trials_NI)-trial_category(2, i))*6;
                    else
                        temp = (trials-trial_category(2, i))*6;
                    end
                end
                plot(hsp{2}, t_trajectory, temp+(z_nose(:, 1)-z_min)*z2y, 'Color', [2/3 2/3 2/3], 'LineWidth', 1);
                plot(hsp{2}, t_trajectory, temp+(z_pawr(:, 1)-z_min)*z2y, 'Color', [1/3 1/3 1/3], 'LineWidth', 1);
                plot(hsp{2}, t_trajectory, temp+(z_pawl(:, 1)-z_min)*z2y, 'Color', [0 0 0], 'LineWidth', 1);
                
                plot(hsp{3}, t_trajectory, temp+(pawl2nose-dz_min)*dz2y, 'Color', [0 0 0], 'LineWidth', 1);
                plot(hsp{3}, t_trajectory, temp+(pawr2nose-dz_min)*dz2y, 'Color', [0.5 0.5 0.5], 'LineWidth', 1);
            end
        end
    end
end

if app.OptSessionCheckBox.Value
    for i = 1:trials
        for ii = 1:4
            if trial_category(1, i)
                LabelledEvents = trials_I(trial_category(2, i)).LabelledEvents;
                bite_timestamps = trials_I(trial_category(2, i)).bite_timestamps;
                if ~isempty(bite_timestamps)
                    bite_amplitudes = trials_I(trial_category(2, i)).bite_amplitudes;
                end
                temp = numel(trials_NI)+trial_category(2, i);
            else
                LabelledEvents = trials_NI(trial_category(2, i)).LabelledEvents;
                bite_timestamps = trials_NI(trial_category(2, i)).bite_timestamps;
                if ~isempty(bite_timestamps)
                    bite_amplitudes = trials_NI(trial_category(2, i)).bite_amplitudes;
                end
                temp = trial_category(2, i);
            end
            
            if ~isempty(LabelledEvents)
                plot_events(hsp{ii}, temp, trials, LabelledEvents);
            end
            
            if ~isempty(bite_timestamps)
                for j = 1:numel(bite_timestamps)
                    line(hsp{ii}, [bite_timestamps(j) bite_timestamps(j)], [(trials-temp)*6+0.5-bite_amplitudes(j)/2 (trials-temp)*6+0.5+bite_amplitudes(j)/2],...
                        'Color', [1 0 1], 'LineStyle', '-', 'LineWidth', 1);
                end
            end
        end
    end
end

if zscore_max > 0 && zscore_min < 0
    zscore2y = 6/(zscore_max-zscore_min);
    zero_position = zscore2y*(-zscore_min);
    for i = 1:trials
        if ~isempty(fpdata_all)
            fpdata = fpdata_all.zsignal_all(:, value(i));
            fpdata_zsignal = [];
            fpdata_t = [];
            for j = 1:nchannel
                fpdata_zsignal(:, j) = fpdata{j}(:, 2);
            end
            fpdata_t = fpdata{1}(:, 1);
            if maxtime(1, i) ~= 0
                fpdata_zsignal(fpdata_t > maxtime(1, i)+4, :) = [];
                fpdata_t(fpdata_t > maxtime(1, i)+4) = [];
            elseif maxtime(2, i) ~= 0
                fpdata_zsignal(fpdata_t > maxtime(2, i)+4, :) = [];
                fpdata_t(fpdata_t > maxtime(2, i)+4) = [];
            end
            maxtime(4, i) = fpdata_t(end);
            
            for j = 1:nchannel
                plot(hsp{4}, fpdata_t, (trials-i)*6+zero_position+fpdata_zsignal(:, j)*zscore2y, 'Color', [1-1/nchannel*j 1-1/nchannel*j 1-1/nchannel*j], 'LineWidth', 1);
            end
            line(hsp{4}, [-4, fpdata_t(end)], [(trials-i)*6+zero_position (trials-i)*6+zero_position],...
                'Color', [0 0 0], 'LineStyle', '--', 'LineWidth', 0.5)
        end
    end
end
for i = 1:4
    if app.OptSessionCheckBox.Value
        laser_timestamps = round(laser_timestamps);
        for j = 1:numel(laser_timestamps)
            line(hsp{i}, [laser_timestamps(j) laser_timestamps(j)], [0 6*trials],...
                'Color', [1 0 0], 'LineStyle', '--', 'LineWidth', 1);
        end
    end
    
    ylim(hsp{i}, [0 trials*6]);
    switch i
        case 1
            if max(maxtime(1, :)) > 0
                xlim(hsp{i}, [0 max(maxtime(1, :))]);
            end
        case {2, 3}
            if max(maxtime(1:3, :), [], 'all') > 0
                xlim(hsp{i}, [0 max(maxtime(1:3, :), [], 'all')]);
            end
        case 4
            line(hsp{i}, [0 0], [0 6*trials],...
                'Color', [0 0 0], 'LineStyle', ':', 'LineWidth', 1)
            xlim(hsp{i}, [-4 max(maxtime([1 2 4], :), [], 'all')]);
    end
    xlabel(hsp{i}, 'Time (s)');
    ylabel(hsp{i}, 'Trial #');
    set(hsp{i}, 'YTick', 3:6:trials*6-3, 'YTickLabel', num2str((trials:-1:1)'), 'TickDir', 'out', 'FontSize', 12, 'ButtonDownFcn', @extract_figure);
end