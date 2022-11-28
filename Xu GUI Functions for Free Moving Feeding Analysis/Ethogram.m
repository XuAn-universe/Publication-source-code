function Ethogram(app, Exp_Path, FrameRate, nmedian)
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

% get photometry data
try
    temp = load([audiolocation '\FPData.mat']);
    fpdata = temp.zsignal_all(:, trial);
    nchannel = numel(fpdata);
    for i = 1:nchannel
        fpdata_zsignal(:, i) = fpdata{i}(:, 2);
        lgdtext{i} = ['Channel ' num2str(i)];
    end
    fpdata_t = fpdata{1}(:, 1);
catch
    fpdata = [];
end

% process trajectories
if app.TrackingDataCheckBox.Value
    label_table = table2array(app.UITable.Data);
else
    helpdlg('''Tracking Data'' needs to be checked');
    return;
end

[~, ~, ~, ~, ~, x_nose, y_nose, z_nose, speed_nose, acceleration_nose, laserstart, laserstop] = trajectory_postprocessing(3, Exp_Path, trial,...
    label_table, nmedian, FrameRate);
[~, ~, ~, ~, ~, x_pawl, y_pawl, z_pawl, speed_pawl, acceleration_pawl, ~, ~] = trajectory_postprocessing(10, Exp_Path, trial, label_table, nmedian, FrameRate);
[~, ~, ~, ~, ~, x_pawr, y_pawr, z_pawr, speed_pawr, acceleration_pawr, ~, ~] = trajectory_postprocessing(16, Exp_Path, trial, label_table, nmedian, FrameRate);

pawl2nose = z_pawl(:, 1)-z_nose(:, 1);
pawr2nose = z_pawr(:, 1)-z_nose(:, 1);
pawl2pawr = z_pawl(:, 1)-z_pawr(:, 1);
t = (1:size(z_nose, 1))'/FrameRate;
time_range = [0 inf];
if app.NoLightButton.Value || ~app.OptSessionCheckBox.Value
    laserstart = [];
    laserstop = [];
else
    laserstart = laserstart/FrameRate;
    laserstop = laserstop/FrameRate;
end

scalefactor = 0.1;
y = [pawl2nose pawr2nose pawl2pawr];
yrange = range(y(:));

% relative distance and velocity
figure
yyaxis left;
hp1 = plot_tj_singletrial(t, y, time_range, laserstart, laserstop,...
    [0 1 0], {'r', 'g', 'b'}, 'Distance (mm)', {'PawL2Nose', 'PawR2Nose', 'PawL2PawR'});

if ~isempty(bite_timestamps)
    add_biteevents(bite_timestamps, [min(y(:)) scalefactor*yrange+min(y(:))]);
end

yyaxis right;
hp2 = plot_tj_singletrial(t, [speed_nose(:, 3) speed_pawl(:, 3) speed_pawr(:, 3)], time_range, laserstart, laserstop,...
    [0 1 0], {'c', 'm', 'y'}, 'Velocity (mm/s)', {'Nose', 'PawL', 'PawR'});

legend off;
legend([hp1 hp2], {'PawL2Nose', 'PawR2Nose', 'PawL2PawR', 'Nose', 'PawL', 'PawR'});

% relative distance and speed
figure
yyaxis left;
hp1 = plot_tj_singletrial(t, y, time_range, laserstart, laserstop,...
    [0 1 0], {'r', 'g', 'b'}, 'Distance (mm)', {'PawL2Nose', 'PawR2Nose', 'PawL2PawR'});

if ~isempty(bite_timestamps)
    add_biteevents(bite_timestamps, [min(y(:)) scalefactor*yrange+min(y(:))]);
end

yyaxis right;
hp2 = plot_tj_singletrial(t, abs([speed_nose(:, 3) speed_pawl(:, 3) speed_pawr(:, 3)]), time_range, laserstart, laserstop,...
    [0 1 0], {'c', 'm', 'y'}, 'Speed (mm/s)', {'Nose', 'PawL', 'PawR'});

legend off;
legend([hp1 hp2], {'PawL2Nose', 'PawR2Nose', 'PawL2PawR', 'Nose', 'PawL', 'PawR'});

% relative distance and acceleration
figure
yyaxis left;
hp1 = plot_tj_singletrial(t, y, time_range, laserstart, laserstop,...
    [0 1 0], {'r', 'g', 'b'}, 'Distance (mm)', {'PawL2Nose', 'PawR2Nose', 'PawL2PawR'});

if ~isempty(bite_timestamps)
    add_biteevents(bite_timestamps, [min(y(:)) scalefactor*yrange+min(y(:))]);
end

yyaxis right;
hp2 = plot_tj_singletrial(t, [acceleration_nose(:, 3) acceleration_pawl(:, 3) acceleration_pawr(:, 3)], time_range, laserstart, laserstop,...
    [0 1 0], {'c', 'm', 'y'}, 'Acceleration (mm/s^{2})', {'Nose', 'PawL', 'PawR'});

legend off;
legend([hp1 hp2], {'PawL2Nose', 'PawR2Nose', 'PawL2PawR', 'Nose', 'PawL', 'PawR'});

% relative distance and absolute acceleration
figure
yyaxis left;
hp1 = plot_tj_singletrial(t, y, time_range, laserstart, laserstop,...
    [0 1 0], {'r', 'g', 'b'}, 'Distance (mm)', {'PawL2Nose', 'PawR2Nose', 'PawL2PawR'});

if ~isempty(bite_timestamps)
    add_biteevents(bite_timestamps, [min(y(:)) scalefactor*yrange+min(y(:))]);
end

yyaxis right;
hp2 = plot_tj_singletrial(t, abs([acceleration_nose(:, 3) acceleration_pawl(:, 3) acceleration_pawr(:, 3)]), time_range, laserstart, laserstop,...
    [0 1 0], {'c', 'm', 'y'}, 'Absolute Acceleration (mm/s^{2})', {'Nose', 'PawL', 'PawR'});

legend off;
legend([hp1 hp2], {'PawL2Nose', 'PawR2Nose', 'PawL2PawR', 'Nose', 'PawL', 'PawR'});

%%
if ~isempty(fpdata)
    y = fpdata_zsignal;
    yrange = range(y(:));
    lcolor = {'r', 'g', 'b'};
    % z position vs gcamp signal
    figure
    yyaxis left;
    hp1 = plot_tj_singletrial(fpdata_t, fpdata_zsignal, time_range, laserstart, laserstop,...
        [0 1 0], lcolor(1:nchannel), 'Z score', lgdtext);
    
    if ~isempty(bite_timestamps)
        add_biteevents(bite_timestamps, [min(y(:)) scalefactor*yrange+min(y(:))]);
    end
    
    yyaxis right;
    hp2 = plot_tj_singletrial(t, [z_nose(:, 1) z_pawl(:, 1) z_pawr(:, 1)], time_range, laserstart, laserstop,...
        [0 1 0], {'c', 'm', 'y'}, 'Z (mm)', {'Nose', 'PawL', 'PawR'});
    
    legend off;
    legend([hp1 hp2], [lgdtext(:); {'Nose'; 'PawL'; 'PawR'}]);
    
    % relative distance vs gcamp signal
    figure
    yyaxis left;
    hp1 = plot_tj_singletrial(fpdata_t, fpdata_zsignal, time_range, laserstart, laserstop,...
        [0 1 0], lcolor(1:nchannel), 'Z score', lgdtext);   
    
    if ~isempty(bite_timestamps)
        add_biteevents(bite_timestamps, [min(y(:)) scalefactor*yrange+min(y(:))]);
    end
    
    yyaxis right;
    hp2 = plot_tj_singletrial(t, [pawl2nose pawr2nose pawl2pawr], time_range, laserstart, laserstop,...
        [0 1 0], {'c', 'm', 'y'}, 'Distance (mm)', {'PawL2Nose', 'PawR2Nose', 'PawL2PawR'});
    
    legend off;
    legend([hp1 hp2], [lgdtext(:); {'PawL2Nose'; 'PawR2Nose'; 'PawL2PawR'}]);
    
    % velocity vs gcamp signal
    figure
    yyaxis left;
    hp1 = plot_tj_singletrial(fpdata_t, fpdata_zsignal, time_range, laserstart, laserstop,...
        [0 1 0], lcolor(1:nchannel), 'Z score', lgdtext);   
    
    if ~isempty(bite_timestamps)
        add_biteevents(bite_timestamps, [min(y(:)) scalefactor*yrange+min(y(:))]);
    end
    
    yyaxis right;
    hp2 = plot_tj_singletrial(t, [speed_nose(:, 3) speed_pawl(:, 3) speed_pawr(:, 3)], time_range, laserstart, laserstop,...
        [0 1 0], {'c', 'm', 'y'}, 'Velocity (mm/s)', {'Nose', 'PawL', 'PawR'});
    
    legend off;
    legend([hp1 hp2], [lgdtext(:); {'Nose'; 'PawL'; 'PawR'}]);
    
    % speed vs gcamp signal
    figure
    yyaxis left;
    hp1 = plot_tj_singletrial(fpdata_t, fpdata_zsignal, time_range, laserstart, laserstop,...
        [0 1 0], lcolor(1:nchannel), 'Z score', lgdtext);   
    
    if ~isempty(bite_timestamps)
        add_biteevents(bite_timestamps, [min(y(:)) scalefactor*yrange+min(y(:))]);
    end
    
    yyaxis right;
    hp2 = plot_tj_singletrial(t, abs([speed_nose(:, 3) speed_pawl(:, 3) speed_pawr(:, 3)]), time_range, laserstart, laserstop,...
        [0 1 0], {'c', 'm', 'y'}, 'Speed (mm/s)', {'Nose', 'PawL', 'PawR'});
    
    legend off;
    legend([hp1 hp2], [lgdtext(:); {'Nose'; 'PawL'; 'PawR'}]);
    
    % relative distance vs z position vs gcamp signal
    figure
    yyaxis left;
    hp1 = plot_tj_singletrial(fpdata_t, fpdata_zsignal, time_range, laserstart, laserstop,...
        [0 1 0], lcolor(1:nchannel), 'Z score', lgdtext);   
    
    if ~isempty(bite_timestamps)
        add_biteevents(bite_timestamps, [min(y(:)) scalefactor*yrange+min(y(:))]);
    end
    
    yyaxis right;
    hp2 = plot_tj_singletrial(t, [pawl2nose z_pawl(:, 1)], time_range, laserstart, laserstop,...
        [0 1 0], {'c', 'm'}, 'mm', {'PawL2Nose', 'PawL'});
    
    legend off;
    legend([hp1 hp2], [lgdtext(:); {'PawL2Nose'; 'PawL'}]);
end
%%
if ~isempty(bite_timestamps)
    frames2remove = (t < bite_timestamps(1) | t > bite_timestamps(end));
    pawl2nose(frames2remove) = [];
    pawr2nose(frames2remove) = [];
    pawl2pawr(frames2remove) = [];
    speed_nose(frames2remove, :) = [];
    speed_pawl(frames2remove, :) = [];
    speed_pawr(frames2remove, :) = [];
    
    %%
    figure;
    plot(pawl2nose, pawr2nose, 'ok');
    line([min(pawl2nose) max(pawl2nose)], [min(pawl2nose) max(pawl2nose)], 'Color', [1 0 0], 'LineStyle', '--', 'LineWidth', 1)
    xlabel('pawl2nose (mm)');
    ylabel('pawr2nose (mm)');
    title('\DeltaZ (first bite to last bite)');
    box off;
    set(gca, 'TickLength', [0 0], 'FontSize', 12);
    
    figure;
    plot(pawl2nose, pawl2pawr, 'ok');
    xlabel('pawl2nose (mm)');
    ylabel('pawl2pawr (mm)');
    title('\DeltaZ (first bite to last bite)');
    box off;
    set(gca, 'TickLength', [0 0], 'FontSize', 12);
    
    figure;
    plot(pawr2nose, pawl2pawr, 'ok');
    xlabel('pawr2nose (mm)');
    ylabel('pawl2pawr (mm)');
    title('\DeltaZ (first bite to last bite)');
    box off;
    set(gca, 'TickLength', [0 0], 'FontSize', 12);
    
    figure;
    plot(pawl2nose, abs(speed_pawl(:, 3)), 'ok');
    xlabel('pawl2nose (mm)');
    ylabel('pawl speed (mm/s)');
    title('(first bite to last bite)');
    box off;
    set(gca, 'TickLength', [0 0], 'FontSize', 12);
    
    figure;
    plot(pawr2nose, abs(speed_pawr(:, 3)), 'ok');
    xlabel('pawr2nose (mm)');
    ylabel('pawr speed (mm/s)');
    title('(first bite to last bite)');
    box off;
    set(gca, 'TickLength', [0 0], 'FontSize', 12);
    %%
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
    plot(hz_filtkern,fft_filtkern(1:ceil(length(fft_filtkern)/2)),'b')
    
    set(gca,'ylim',[-.1 1.1],'xlim',[0 nyquist])
    legend({'ideal';'best fit'})
    
    freqsidx = dsearchn(hz_filtkern',ffrequencies'*nyquist);
    title([ 'SSE: ' num2str(sum( (idealresponse-fft_filtkern(freqsidx)).^2 )) ]);
    
    z_nose(frames2remove, :) = [];
    z_nose = fillmissing(z_nose, 'movmedian', nmedian);
    z_pawl(frames2remove, :) = [];
    z_pawl = fillmissing(z_pawl, 'movmedian', nmedian);
    z_pawr(frames2remove, :) = [];
    z_pawr = fillmissing(z_pawr, 'movmedian', nmedian);
    t(frames2remove, :) = [];
    filter_znose = filtfilt(filterweights, 1, z_nose(:, 1)); % apply filter to the data
    filter_zpawl = filtfilt(filterweights, 1, z_pawl(:, 1)); % apply filter to the data
    filter_zpawr = filtfilt(filterweights, 1, z_pawr(:, 1)); % apply filter to the data
    scalefactor = 0.1;
    y = [filter_znose; filter_zpawl; filter_zpawr];
    yrange = range(y(:));
    figure;
    hp(1) = plot(t, filter_znose, '-b');
    hold on;
    hp(2) = plot(t, filter_zpawl, '-r');
    hp(3) = plot(t, filter_zpawr, '-g');
    if ~isempty(bite_timestamps)
        add_biteevents(bite_timestamps, [min(y(:)) scalefactor*yrange+min(y(:))]);
    end
    plot([bite_timestamps(1) bite_timestamps(end)], [0 0], '-k');
    xlabel('Time (s)');
    ylabel('Z (mm)');
    title('Z (first bite to last bite)');
    box off;
    if ~isempty(laserstart)
        for j = 1:numel(laserstart)
            %         patch([laserstart(j) laserstart(j) laserstop(j) laserstop(j)],...
            %             [min(y(:))-scalefactor*yrange max(y(:))+scalefactor*yrange max(y(:))+scalefactor*yrange min(y(:))-scalefactor*yrange],...
            %             pcolor, 'FaceAlpha', 1, 'EdgeColor', 'none');
            
            patch([laserstart(j) laserstart(j) laserstop(j) laserstop(j)],...
                [max(y(:))+scalefactor*yrange max(y(:))+2*scalefactor*yrange max(y(:))+2*scalefactor*yrange max(y(:))+scalefactor*yrange],...
                [0 1 0], 'FaceAlpha', 1, 'EdgeColor', 'none');
        end
    end
    ylim([min(y(:)) max(y(:))+2*scalefactor*yrange]);
    set(gca, 'TickLength', [0 0], 'FontSize', 12);
    legend(hp, {'Nose', 'Left Paw', 'Right Paw'}, 'Location', 'northeast', 'FontSize', 12)
end



function add_biteevents(bite_timestamps, y)
if nargin == 1
    y = ylim;
end
for i = 1:numel(bite_timestamps)
    line([bite_timestamps(i) bite_timestamps(i)], y,...
        'Color', [0 0 0], 'LineStyle', '-', 'LineWidth', 1);
end