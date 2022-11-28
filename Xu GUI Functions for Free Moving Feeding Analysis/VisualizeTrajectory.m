function VisualizeTrajectory(app, Exp_Path, FrameRate, nmedian, fpdata)
filtering4phase = 0; % whether to filter z trajectory to compute phase
trial = app.TrialsListBox.Value;
if numel(trial) ~= 1
    helpdlg('You can only choose one trial');
    return;
end
load([Exp_Path '\Analysis_Session.mat'], 'Video_annotation');
time_ready2bite = Video_annotation(trial).time_ready2bite;
time_feeding_end = Video_annotation(trial).time_feeding_end;

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
% get photometry data
if ~isempty(fpdata)
    nchannel = numel(fpdata);
    for i = 1:nchannel
        fpdata_zsignal(:, i) = fpdata{i}(:, 2);
        lgdtext{i} = ['Channel ' num2str(i)];
    end
    fpdata_t = fpdata{1}(:, 1);
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

t = (1:numel(y_PG1))'/FrameRate;
time_range = [app.timerangesEditField.Value app.toEditField.Value];
if app.NoLightButton.Value || ~app.OptSessionCheckBox.Value
    laserstart = [];
    laserstop = [];
else
    laserstart = laserstart/FrameRate;
    laserstop = laserstop/FrameRate;
end

% raw X, Y, Z trajectories
figure;
subplot(3, 3, 1);
plot_tj_singletrial(t, x(:, 1), time_range, laserstart, laserstop, [0 1 0], {'k'}, 'X (mm)', {'X'});

subplot(3, 3, 2);
plot_tj_singletrial(t, [y(:, 1) y_PG1 y_PG3], time_range, laserstart, laserstop,...
    [0 1 0], {'k', 'r', 'b'}, 'Y (mm)', {'Y', 'Y PG1', 'Y PG3'});

subplot(3, 3, 3);
plot_tj_singletrial(t, [z(:, 1) z_PG1 z_PG2 z_PG3], time_range, laserstart, laserstop, [0 1 0], {'k', 'r', 'm', 'b'}, 'Z (mm)',...
    {'Z', 'Z PG1', 'Z PG2', 'Z PG3'});

% raw and filtered X, Y, Z trajectories
subplot(3, 3, 4);
plot_tj_singletrial(t, x, time_range, laserstart, laserstop, [0 1 0], {'r', 'k'}, 'X (mm)', {'X', 'X filtered'});

subplot(3, 3, 5);
plot_tj_singletrial(t, y, time_range, laserstart, laserstop, [0 1 0], {'r', 'k'}, 'Y (mm)', {'Y', 'Y filtered'});

subplot(3, 3, 6);
plot_tj_singletrial(t, z, time_range, laserstart, laserstop, [0 1 0], {'r', 'k'}, 'Z (mm)', {'Z', 'Z filtered'});

% speed
subplot(3, 3, 7);
plot_tj_singletrial(t, speed, time_range, laserstart, laserstop, [0.9 0.9 0.9], {'r', 'g', 'b', 'y', 'm', 'c', 'k'}, 'Velocity or Speed (mm/s)',...
    {'X', 'Y', 'Z', 'XY', 'XZ', 'YZ', 'XYZ'});

% acceleration
subplot(3, 3, 8);
plot_tj_singletrial(t, acceleration, time_range, laserstart, laserstop, [0.9 0.9 0.9], {'r', 'g', 'b', 'y', 'm', 'c', 'k'},...
    'Acceleration (mm/s^{2})', {'X', 'Y', 'Z', 'XY', 'XZ', 'YZ', 'XYZ'});

% raw X, Y, Z trajectories with bite events
figure;
subplot(1, 3, 1);
plot_tj_singletrial(t, x(:, 1), time_range, laserstart, laserstop, [0 1 0], {'k'}, 'X (mm)', {'X'});
if ~isempty(bite_timestamps)
    plot_biteevents(bite_timestamps, bite_amplitudes);
end

subplot(1, 3, 2);
plot_tj_singletrial(t, y(:, 1), time_range, laserstart, laserstop, [0 1 0], {'k'}, 'Y (mm)', {'Y'});
if ~isempty(bite_timestamps)
    plot_biteevents(bite_timestamps, bite_amplitudes);
end

subplot(1, 3, 3);
plot_tj_singletrial(t, z(:, 1), time_range, laserstart, laserstop, [0 1 0], {'k'}, 'Z (mm)', {'Z'});
if ~isempty(bite_timestamps)
    plot_biteevents(bite_timestamps, bite_amplitudes);
end

if numel(bite_timestamps) >= 2 || (~isempty(time_ready2bite) && ~isempty(time_feeding_end) && isempty(bite_timestamps))
    if isempty(bite_timestamps)
        bite_timestamps(1) = time_ready2bite;
        bite_timestamps(2) = time_feeding_end;
    end
    t_bite = t(t >= bite_timestamps(1) & t <= bite_timestamps(end));
    z_bite = z(t >= bite_timestamps(1) & t <= bite_timestamps(end), 1);
    try
        % Fourier transform of Z
        [f, MX] = FFT_trajectory(fillmissing(z_bite, 'movmedian', nmedian), FrameRate/2);
        figure;
        plot(f, sqrt(MX), '-k');
        xlim([0 FrameRate/2]);
        xlabel('Frequency (Hz)');
        ylabel('Amplitude (mm)');
        title('Z (first bite to last bite)');
        set(gca, 'TickLength', [0 0], 'FontSize', 12);
        box off;
        
        if filtering4phase
            % 1-0. Initiate the parameter
            nyquist = FrameRate / 2; % theoretical limit
            pband1 = 0.15;
            pband2 = 2;
            filterlen = round(1/pband1*FrameRate*1);
            
            % 1-1. Creat a Delat Filter (from DeltaFilter4sec)
            transition_width = 0.1; % usually 0.1~0.25, the sharpeness of the filter
            
            ffrequencies   = [ 0 (1-transition_width)*pband1 pband1 pband2 (1+transition_width)*pband2 nyquist ]/nyquist;
            idealresponse  = [ 0 0 1 1 0 0 ]; % Band-pass filter
            filterweights = firls(filterlen,ffrequencies,idealresponse);
            
            % 1-2. Check the Filter Quality
            filterweightsW = zscore(filterweights); % with z-score
            figure; plot(ffrequencies*nyquist,idealresponse,'r');hold on
            
            fft_filtkern  = abs(fft(filterweightsW));
            fft_filtkern  = fft_filtkern./max(fft_filtkern); % normalized to 1.0 for visual comparison ease
            
            hz_filtkern =linspace(0, nyquist, floor(length(fft_filtkern)/2+1));
            plot(hz_filtkern,fft_filtkern(1:ceil(length(fft_filtkern)/2)),'b');
            
            set(gca,'ylim',[-.1 1.1],'xlim',[0 nyquist]);
            set(gca, 'TickLength', [0 0], 'FontSize', 12);
            legend({'ideal';'best fit'});
            
            freqsidx = dsearchn(hz_filtkern',ffrequencies'*nyquist);
            title([ 'SSE: ' num2str(sum((idealresponse-fft_filtkern(freqsidx)).^2 )) ]);
            
            filter_data = filtfilt(filterweights, 1, fillmissing(z_bite, 'movmedian', nmedian)); % apply filter to the data
            figure;
            hp(1) = plot(t_bite, filter_data, '-b');
            hold on;
            if ~isempty(bite_amplitudes)
                plot_biteevents(bite_timestamps, bite_amplitudes);
            end
            hp(2) = plot(t_bite, z_bite, '-k');
            plot([bite_timestamps(1) bite_timestamps(end)], [0 0], '-k');
            xlabel('Time (s)');
            ylabel('Z (mm)');
            title('Z (first bite to last bite)');
            set(gca, 'TickLength', [0 0], 'FontSize', 12);
            box off;
            yl = ylim;
            if ~isempty(laserstart)
                for j = 1:numel(laserstart)
                    line([laserstart(j) laserstart(j)], yl, 'Color', [0 1 0], 'LineStyle', '--', 'LineWidth', 1);
                    line([laserstop(j) laserstop(j)], yl, 'Color', [0 1 0], 'LineStyle', '--', 'LineWidth', 1);
                end
            end
            legend(hp, {'band pass filtered', 'original'}, 'Location', 'northeast', 'FontSize', 12);
        end
        
        if ~isempty(bite_amplitudes)
            if filtering4phase
                hilbert_data = hilbert(filter_data);
            else
                hilbert_data = hilbert(fillmissing(z_bite, 'movmedian', nmedian));
            end
            phase_data = angle(hilbert_data)/pi*180; % this angle is the cosine angle
            figure;
            plot(t_bite, phase_data, '-k');
            hold on;
            phase_bite = plot_biteevents(bite_timestamps, bite_amplitudes);
            xlabel('Time (s)');
            ylabel('Phase');
            title('Z (first bite to last bite)');
            set(gca, 'TickLength', [0 0], 'FontSize', 12);
            box off;
            yl = ylim;
            if ~isempty(laserstart)
                for j = 1:numel(laserstart)
                    line([laserstart(j) laserstart(j)], yl, 'Color', [0 1 0], 'LineStyle', '--', 'LineWidth', 1);
                    line([laserstop(j) laserstop(j)], yl, 'Color', [0 1 0], 'LineStyle', '--', 'LineWidth', 1);
                end
            end
            figure;
            histogram(phase_bite, -180:180/12:180, 'FaceColor', [0 0 0], 'FaceAlpha', 1, 'Normalization', 'probability');
            xlim([-180 180]);
            set(gca, 'TickLength', [0 0], 'FontSize', 12);
            box off;
            xlabel('Phase');
            ylabel('probability');
        end
    end
    
    if ~isempty(fpdata)
        try
            % 1-0. Initiate the parameter
            SampleRate = 1/mean(diff(fpdata_t));
            nyquist = SampleRate / 2; % theoretical limit
            filterlen = round(1/pband1*SampleRate*1);
            
            % 1-1. Creat a Delat Filter (from DeltaFilter4sec)
            transition_width = 0.1; % usually 0.1~0.25, the sharpeness of the filter
            
            ffrequencies   = [ 0 (1-transition_width)*pband1 pband1 pband2 (1+transition_width)*pband2 nyquist ]/nyquist;
            idealresponse  = [ 0 0 1 1 0 0 ]; % Band-pass filter
            filterweights = firls(filterlen,ffrequencies,idealresponse);
            
            % 1-2. Check the Filter Quality
            filterweightsW = zscore(filterweights); % with z-score
            figure; plot(ffrequencies*nyquist,idealresponse,'r');hold on
            
            fft_filtkern  = abs(fft(filterweightsW));
            fft_filtkern  = fft_filtkern./max(fft_filtkern); % normalized to 1.0 for visual comparison ease
            
            hz_filtkern =linspace(0, nyquist, floor(length(fft_filtkern)/2+1));
            plot(hz_filtkern,fft_filtkern(1:ceil(length(fft_filtkern)/2)),'b');
            
            set(gca,'ylim',[-.1 1.1],'xlim',[0 nyquist]);
            set(gca, 'TickLength', [0 0], 'FontSize', 12);
            legend({'ideal';'best fit'});
            
            freqsidx = dsearchn(hz_filtkern',ffrequencies'*nyquist);
            title([ 'SSE: ' num2str(sum((idealresponse-fft_filtkern(freqsidx)).^2 )) ]);
            
            figure;
            for i = 1:nchannel
                filter_fpdata(:, i) = filtfilt(filterweights, 1, fpdata_zsignal(fpdata_t >= bite_timestamps(1) & fpdata_t <= bite_timestamps(end), i)); % apply filter to the data
                subplot(1, nchannel, i);
                hp(1) = plot(fpdata_t(fpdata_t >= bite_timestamps(1) & fpdata_t <= bite_timestamps(end)), fpdata_zsignal(fpdata_t >= bite_timestamps(1) & fpdata_t <= bite_timestamps(end), i), '-k');
                hold on;
                if ~isempty(bite_amplitudes)
                    plot_biteevents(bite_timestamps, bite_amplitudes);
                end
                hp(2) = plot(fpdata_t(fpdata_t >= bite_timestamps(1) & fpdata_t <= bite_timestamps(end)), filter_fpdata(:, i), '-b');
                plot([bite_timestamps(1) bite_timestamps(end)], [0 0], '-k');
                xlabel('Time (s)');
                ylabel('Z score');
                title(['Z score (first bite to last bite): Channel ' num2str(i)]);
                set(gca, 'TickLength', [0 0], 'FontSize', 12);
                box off;
                legend(hp, {'original', 'band pass filtered'}, 'Location', 'northeast', 'FontSize', 12);
            end
            
            for i = 1:nchannel
                figure;
                hold on;
                yyaxis left;
                hp(1) = plot(fpdata_t(fpdata_t >= bite_timestamps(1) & fpdata_t <= bite_timestamps(end)), filter_fpdata(:, i), '-g');
                ylabel('Z score');
                set(gca, 'YColor', [0 0.5 0]);
                
                yyaxis right;
                hp(2) = plot(t_bite, filter_data, '-k');
                plot([bite_timestamps(1) bite_timestamps(end)], [0 0], '-k');
                xlabel('Time (s)');
                ylabel('Z (mm)');
                title(['Z score vs Z position (first bite to last bite): Channel ' num2str(i)]);
                set(gca, 'TickLength', [0 0], 'FontSize', 12, 'YColor', [0 0 0]);
                box off;
                legend(hp, {'Z score', 'Z position'}, 'Location', 'northeast', 'FontSize', 12);
            end
        end
        
        win = 5;
        z_smooth = smoothdata(z(:, 1), 'gaussian', win, 'omitnan');
        z_bite_smooth = z_smooth(t >= bite_timestamps(1) & t <= bite_timestamps(end));
        %         figure;
        %         plot(t, z(:, 1), '-k');
        %         hold on;
        %         plot(t, z_smooth, '-r');
        
        zspeed = [NaN; diff(z_smooth)];
        zspeed_bite = zspeed(t >= bite_timestamps(1) & t <= bite_timestamps(end));
        if any(isnan(zspeed_bite))
            errordlg('There is NaN in the velocity varible.', 'Error');
            return;
        end
        [risingindex, fallingindex, rise2fall, fall2rise] = ZeroCrossingDetection(zspeed_bite);
        
        figure;
        hold on;
        yyaxis left;
        hp(1) = plot(t_bite, z_bite_smooth, '-k');
        xlabel('Time (s)');
        ylabel('Z (mm)');
        set(gca, 'YColor', [0 0 0]);
        
        yyaxis right;
        hp(2) = plot(t_bite, zspeed_bite, '-b');
        plot(t_bite(risingindex), zspeed_bite(risingindex), 'og');
        plot(t_bite(fallingindex), zspeed_bite(fallingindex), 'or');
        ylabel('Velocity (mm/s)');
        line([bite_timestamps(1) bite_timestamps(end)], [0 0], 'Color', [0 0 0],  'LineStyle', '-', 'LineWidth', 1);
        xlim([bite_timestamps(1) bite_timestamps(end)]);
        set(gca, 'YColor', [0 0 1]);
        set(gca, 'TickLength', [0 0], 'FontSize', 12);
        title('Zero crossings of velocity');
        legend(hp, {'Z position', 'Z velocity'});
        
        riseheight = z_bite_smooth(rise2fall(:, 2))-z_bite_smooth(rise2fall(:, 1));
        fallheight = z_bite_smooth(fall2rise(:, 1))-z_bite_smooth(fall2rise(:, 2));
        threshold = 0.1; % 5%
        [rise_height_top, rise_index_top] = maxk(riseheight, round(threshold*numel(riseheight)));
        [fall_height_top, fall_index_top] = maxk(fallheight, round(threshold*numel(fallheight)));
        figure;
        subplot(2, 2, 1);
        histogram(riseheight, 'FaceColor', [0 0 0], 'FaceAlpha', 1, 'Normalization', 'probability');
        yl = ylim;
        line([min(rise_height_top) min(rise_height_top)], yl, 'Color', [1 0 0], 'LineStyle', '-', 'LineWidth', 1)
        set(gca, 'TickLength', [0 0], 'FontSize', 12, 'ButtonDownFcn', @extract_figure);
        box off;
        xlabel('\DeltaZ (mm)');
        ylabel('probability');
        title('lift height distribution');
        
        subplot(2, 2, 2);
        histogram(fallheight, 'FaceColor', [0 0 0], 'FaceAlpha', 1, 'Normalization', 'probability');
        yl = ylim;
        line([min(fall_height_top) min(fall_height_top)], yl, 'Color', [1 0 0], 'LineStyle', '-', 'LineWidth', 1)
        set(gca, 'TickLength', [0 0], 'FontSize', 12, 'ButtonDownFcn', @extract_figure);
        box off;
        xlabel('\DeltaZ (mm)');
        ylabel('probability');
        title('fall height distribution');
        
        subplot(2, 2, 3);
        hold on;
        plot(t_bite, z_bite_smooth, '-k');
        if ~isempty(bite_amplitudes)
            plot_biteevents(bite_timestamps, bite_amplitudes);
        end
        plot(t_bite(rise2fall(rise_index_top, 1)), z_bite_smooth(rise2fall(rise_index_top, 1)), 'og');
        plot(t_bite(rise2fall(rise_index_top, 2)), z_bite_smooth(rise2fall(rise_index_top, 2)), 'ob');
        xlabel('Time (s)');
        ylabel('Z (mm)');
        title(['top ' num2str(threshold*100) '% lift height']);
        set(gca, 'TickLength', [0 0], 'FontSize', 12, 'ButtonDownFcn', @extract_figure);
        
        subplot(2, 2, 4);
        hold on;
        plot(t_bite, z_bite_smooth, '-k');
        if ~isempty(bite_amplitudes)
            plot_biteevents(bite_timestamps, bite_amplitudes);
        end
        plot(t_bite(fall2rise(fall_index_top, 1)), z_bite_smooth(fall2rise(fall_index_top, 1)), 'og');
        plot(t_bite(fall2rise(fall_index_top, 2)), z_bite_smooth(fall2rise(fall_index_top, 2)), 'ob');
        xlabel('Time (s)');
        ylabel('Z (mm)');
        title(['top ' num2str(threshold*100) '% fall height']);
        set(gca, 'TickLength', [0 0], 'FontSize', 12, 'ButtonDownFcn', @extract_figure);
        
        pre = app.presEditField.Value;
        post = app.postsEditField.Value;
        [t_lift, zscore_lift, zscore_lift_random] = AlignSignal2Event(fpdata_t, fpdata_zsignal, t_bite(rise2fall(rise_index_top, 1)), pre, post);
        if ~isempty(bite_amplitudes)
            nrise = numel(rise_index_top);
            firstbite = zeros(nrise, 1);
            bitealign2rise = repmat(bite_timestamps, nrise, 1)-t_bite(rise2fall(rise_index_top, 1));
            for i = 1:nrise
                temp = bitealign2rise(i, :);
                firstbite(i) = temp(find(temp >= 0, 1));
            end
            [t_firstbite, zscore_firstbite, zscore_firstbite_random] = AlignSignal2Event(fpdata_t, fpdata_zsignal, t_bite(rise2fall(rise_index_top, 1))+firstbite, pre, post);
        end
        
        [t_fall, zscore_fall, zscore_fall_random] = AlignSignal2Event(fpdata_t, fpdata_zsignal, t_bite(fall2rise(fall_index_top, 1)), pre, post);
        if ~isempty(bite_amplitudes)
            nfall = numel(fall_index_top);
            lastbite = zeros(nfall, 1);
            bitealign2fall = repmat(bite_timestamps, nfall, 1)-t_bite(fall2rise(fall_index_top, 1));
            for i = 1:nfall
                temp = bitealign2fall(i, :);
                lastbite(i) = temp(find(temp >= 0, 1));
            end
            t_lastbite = t_bite(fall2rise(fall_index_top, 1))+lastbite;
            tlimit = 0.1;
            t_lastbite(lastbite > tlimit) = [];
            [t_lastbite, zscore_lastbite, zscore_lastbite_random] = AlignSignal2Event(fpdata_t, fpdata_zsignal, t_lastbite, pre, post);
        end
        
        for i = 1:nchannel
            figure;
            if ~isempty(bite_amplitudes)
                subplot(2, 2, 1);
                yyaxis left;
                plot_gcamp_align2event(t_lift, zscore_lift_random(:, :, i), zscore_lift(:, :, i), ['aligned to lift: Channel ' num2str(i)]);
                set(gca, 'YColor', [0 0.5 0]);
                
                yyaxis right;
                h1 = cdfplot(firstbite);
                set(h1, 'Color', [1 0 0], 'LineWidth', 1.5);
                set(gca, 'GridLineStyle', ':', 'YColor', [1 0 0]);
                xlabel('Time (s)');
                ylabel('Cumulative Probability');
                title(['aligned to lift: Channel ' num2str(i)]);
                
                subplot(2, 2, 2);
                plot_gcamp_align2event(t_firstbite, zscore_firstbite_random(:, :, i), zscore_firstbite(:, :, i), ['aligned to first bite after lift: Channel ' num2str(i)]);
                set(gca, 'ButtonDownFcn', @extract_figure);
                
                subplot(2, 2, 3);
                yyaxis left;
                plot_gcamp_align2event(t_fall, zscore_fall_random(:, :, i), zscore_fall(:, :, i), ['aligned to fall: Channel ' num2str(i)]);
                set(gca, 'YColor', [0 0.5 0]);
                
                yyaxis right;
                h1 = cdfplot(lastbite);
                set(h1, 'Color', [1 0 0], 'LineWidth', 1.5);
                line([tlimit tlimit], [0 1], 'Color', [1 0 0], 'LineStyle', '--', 'LineWidth', 1)
                set(gca, 'GridLineStyle', ':', 'YColor', [1 0 0]);
                xlabel('Time (s)');
                ylabel('Cumulative Probability');
                title(['aligned to fall: Channel ' num2str(i)]);
                
                subplot(2, 2, 4);
                plot_gcamp_align2event(t_lastbite, zscore_lastbite_random(:, :, i), zscore_lastbite(:, :, i), ['aligned to first bite after fall within ' num2str(tlimit) ' s : Channel ' num2str(i)]);
                set(gca, 'ButtonDownFcn', @extract_figure);
            else
                subplot(1, 2, 1);
                plot_gcamp_align2event(t_lift, zscore_lift_random(:, :, i), zscore_lift(:, :, i), ['aligned to lift: Channel ' num2str(i)]);
                set(gca, 'ButtonDownFcn', @extract_figure);
                
                subplot(1, 2, 2);
                plot_gcamp_align2event(t_fall, zscore_fall_random(:, :, i), zscore_fall(:, :, i), ['aligned to fall: Channel ' num2str(i)]);
                set(gca, 'ButtonDownFcn', @extract_figure);
            end
        end
    end
end

% trajectory in 3D
figure;
x = x(t >= time_range(1) & t <= time_range(2), :);
y = y(t >= time_range(1) & t <= time_range(2), :);
z = z(t >= time_range(1) & t <= time_range(2), :);
t = t(t >= time_range(1) & t <= time_range(2));

startID = find(~isnan(x(:, 1)) & ~isnan(y(:, 1)) & ~isnan(z(:, 1)), 1);
endID = find(~isnan(x(:, 1)) & ~isnan(y(:, 1)) & ~isnan(z(:, 1)), 1, 'last');

plot3(x(:, 1), y(:, 1), z(:, 1), 'Color', [0 0 0]);
hold on;
h1 = plot3(x(startID, 1), y(startID, 1), z(startID, 1), 'o', 'Color', [0 0 1], 'MarkerFaceColor', [0 0 1]);
h2 = plot3(x(endID, 1), y(endID, 1), z(endID, 1), 'o', 'Color', [1 0 0], 'MarkerFaceColor', [1 0 0]);
if ~isempty(laserstart)
    for i = 1:numel(laserstart)
        try
            plot3(x(t >= laserstart(i) & t <= laserstop(i), 1), y(t >= laserstart(i) & t <= laserstop(i), 1),...
                z(t >= laserstart(i) & t <= laserstop(i), 1), 'Color', [0 1 0]);
        end
    end
end
xlabel('X (mm)');
ylabel('Y (mm)');
zlabel('Z (mm)');
legend([h1 h2], {'Start', 'End'}, 'Location', 'northeast', 'FontSize', 12);
legend('boxoff');
set(gca, 'TickLength', [0 0], 'FontSize', 12);