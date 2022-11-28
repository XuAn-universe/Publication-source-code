function GCaMP2SingleMovement(app, Exp_Path, FrameRate, nmedian)
Exp_Paths{1} = Exp_Path;

value = app.ListBox.Value;
if ~isempty(value)
    for i = 1:numel(value)
        Exp_Paths{end+1} = app.ListBox.Items{value(i)};
    end
end

if app.TrackingDataCheckBox.Value
    label_table = table2array(app.UITable.Data);
else
    helpdlg('''Tracking Data'' needs to be checked');
    return;
end
pointID = app.DropDown.Value;

window_smooth = 5;
pthreshold = 0.1;
t_limit = 0.1;
pre = app.preEditField.Value;
post = app.postEditField.Value;

ccstep = 4; % number of data point
shifttime = 30; % half of shift time
cc_all = [];
cc2nose_all = [];

riseheight_all = [];
fallheight_all = [];
zscore_lift_all = [];
zscore_lift_random_all = [];
zscore_firstbite_all = [];
zscore_firstbite_random_all = [];
zscore_fall_all = [];
zscore_fall_random_all = [];
zscore_lastbite_all = [];
zscore_lastbite_random_all = [];
firstbite_all = [];
lastbite_all = [];

for i = 1:numel(Exp_Paths)
    load([Exp_Paths{i} '\Analysis_Session.mat'], 'Video_annotation');
    
    % get bite events
    try
        audiolocation = Exp_Paths{i}(1:end-7);
        temp = load([audiolocation '\Detected_Bite_Events.mat']);
        Bite_events = temp.Audio_analysis;
    catch
        Bite_events = [];
    end
    
    for j = 1:numel(Video_annotation)
        if ~Video_annotation(j).Disgard
            [~, ~, ~, ~, ~, ~, ~, z, ~, ~, ~, ~] = trajectory_postprocessing(pointID, Exp_Paths{i}, j, label_table, nmedian, FrameRate);
            t = (1:size(z, 1))'/FrameRate;
            
            if ~isempty(Bite_events)
                bite_timestamps = Bite_events(j).time_bites;
            else
                if ~isempty(Video_annotation(j).time_ready2bite) && ~isempty(Video_annotation(j).time_feeding_end)
                    bite_timestamps(1) = Video_annotation(j).time_ready2bite;
                    bite_timestamps(2) = Video_annotation(j).time_feeding_end;
                end
            end
            
            if numel(bite_timestamps) >= 2
                z_smooth = smoothdata(z(:, 1), 'gaussian', window_smooth, 'omitnan');
                z_bite_smooth = z_smooth(t >= bite_timestamps(1) & t <= bite_timestamps(end));
                
                zspeed = [NaN; diff(z_smooth)];
                zspeed_bite = zspeed(t >= bite_timestamps(1) & t <= bite_timestamps(end));
                if any(isnan(zspeed_bite))
                    disp(['NaN is found in the velocity variable of trial ' num2str(j) ' of ' Exp_Paths{i}]);
                    continue;
                end
                [~, ~, rise2fall, fall2rise] = ZeroCrossingDetection(zspeed_bite);
                riseheight = z_bite_smooth(rise2fall(:, 2))-z_bite_smooth(rise2fall(:, 1));
                fallheight = z_bite_smooth(fall2rise(:, 1))-z_bite_smooth(fall2rise(:, 2));
                riseheight_all = [riseheight_all; riseheight];
                fallheight_all = [fallheight_all; fallheight];
            end
        end
    end
end
[rise_height_top, ~] = maxk(riseheight_all, round(pthreshold*numel(riseheight_all)));
[fall_height_top, ~] = maxk(fallheight_all, round(pthreshold*numel(fallheight_all)));
thr_rise = min(rise_height_top);
thr_fall = min(fall_height_top);

figure;
subplot(1, 2, 1);
histogram(riseheight_all, 'FaceColor', [0 0 0], 'FaceAlpha', 1, 'Normalization', 'probability');
yl = ylim;
line([thr_rise thr_rise], yl, 'Color', [1 0 0], 'LineStyle', '-', 'LineWidth', 1)
set(gca, 'TickLength', [0 0], 'FontSize', 12, 'ButtonDownFcn', @extract_figure);
box off;
xlabel('\DeltaZ (mm)');
ylabel('probability');
title('lift height distribution');

subplot(1, 2, 2);
histogram(fallheight_all, 'FaceColor', [0 0 0], 'FaceAlpha', 1, 'Normalization', 'probability');
yl = ylim;
line([thr_fall thr_fall], yl, 'Color', [1 0 0], 'LineStyle', '-', 'LineWidth', 1)
set(gca, 'TickLength', [0 0], 'FontSize', 12, 'ButtonDownFcn', @extract_figure);
box off;
xlabel('\DeltaZ (mm)');
ylabel('probability');
title('fall height distribution');

for i = 1:numel(Exp_Paths)
    load([Exp_Paths{i} '\Analysis_Session.mat'], 'Video_annotation');
    
    % get bite events
    try
        audiolocation = Exp_Paths{i}(1:end-7);
        temp = load([audiolocation '\Detected_Bite_Events.mat']);
        Bite_events = temp.Audio_analysis;
    catch
        Bite_events = [];
    end
    
    % get photometry data
    try
        fpdata_all = load([audiolocation '\FPData.mat']);
        nchannel = size(fpdata_all.zsignal_all, 1);
        for j = 1:nchannel
            lgdtext{j} = ['Channel ' num2str(j)];
        end
    catch
        errordlg('Fiber photometry data is missing!', 'Error');
        return;
    end
    
    for j = 1:numel(Video_annotation)
        if ~Video_annotation(j).Disgard
            [~, ~, ~, ~, ~, ~, ~, z, ~, ~, ~, ~] = trajectory_postprocessing(pointID, Exp_Paths{i}, j, label_table, nmedian, FrameRate);
            [~, ~, ~, ~, ~, ~, ~, z_nose, ~, ~, ~, ~] = trajectory_postprocessing(3, Exp_Paths{i}, j, label_table, nmedian, FrameRate);
            t = (1:size(z, 1))'/FrameRate;
            
            if ~isempty(Bite_events)
                bite_timestamps = Bite_events(j).time_bites;
            else
                if ~isempty(Video_annotation(j).time_ready2bite) && ~isempty(Video_annotation(j).time_feeding_end)
                    bite_timestamps(1) = Video_annotation(j).time_ready2bite;
                    bite_timestamps(2) = Video_annotation(j).time_feeding_end;
                end
            end
            
            fpdata = fpdata_all.zsignal_all(:, j);
            fpdata_zsignal = [];
            fpdata_t = [];
            for k = 1:nchannel
                fpdata_zsignal(:, k) = fpdata{k}(:, 2);
            end
            fpdata_t = fpdata{1}(:, 1);
            SampleRate = mean(diff(fpdata_t));
            
            if numel(bite_timestamps) >= 2
                t_bite = t(t >= bite_timestamps(1) & t <= bite_timestamps(end));
                z_smooth = smoothdata(z(:, 1), 'gaussian', window_smooth, 'omitnan');
                z_bite_smooth = z_smooth(t >= bite_timestamps(1) & t <= bite_timestamps(end));
                
                z2nose = z(:, 1)-z_nose(:, 1);
                fpdata_zsignal_bite = fpdata_zsignal(fpdata_t >= bite_timestamps(1) & fpdata_t <= bite_timestamps(end), :);
                fpdata_t_bite = fpdata_t(fpdata_t >= bite_timestamps(1) & fpdata_t <= bite_timestamps(end));
                temp = repmat(t_bite', numel(fpdata_t_bite), 1)-fpdata_t_bite;
                [~, index_match] = min(abs(temp), [], 2);
                % compute correlation coefficients
                cc = zeros(2*shifttime+1, 1, nchannel);
                cc2nose = zeros(2*shifttime+1, 1, nchannel);
                for k = -shifttime:1:shifttime
                    z_shift = circshift(z(:, 1), k*ccstep);
                    z_bite_shift = z_shift(t >= bite_timestamps(1) & t <= bite_timestamps(end));
                    z_bite_match = z_bite_shift(index_match);
                    z2nose_shift = circshift(z2nose, k*ccstep);
                    z2nose_bite_shift = z2nose_shift(t >= bite_timestamps(1) & t <= bite_timestamps(end));
                    z2nose_bite_match = z2nose_bite_shift(index_match);
                    for ii = 1:nchannel
                        R = corrcoef(z_bite_match, fpdata_zsignal_bite(:, ii), 'Rows', 'pairwise');
                        cc(k+shifttime+1, 1, ii) = R(2);
                        
                        R = corrcoef(z2nose_bite_match, fpdata_zsignal_bite(:, ii), 'Rows', 'pairwise');
                        cc2nose(k+shifttime+1, 1, ii) = R(2);
                    end
                end
                cc_all = [cc_all cc];
                cc2nose_all = [cc2nose_all cc2nose];
                
                zspeed = [NaN; diff(z_smooth)];
                zspeed_bite = zspeed(t >= bite_timestamps(1) & t <= bite_timestamps(end));
                if any(isnan(zspeed_bite))
                    continue;
                end
                [~, ~, rise2fall, fall2rise] = ZeroCrossingDetection(zspeed_bite);
                
                riseheight = z_bite_smooth(rise2fall(:, 2))-z_bite_smooth(rise2fall(:, 1));
                fallheight = z_bite_smooth(fall2rise(:, 1))-z_bite_smooth(fall2rise(:, 2));
                t_rise = t_bite(rise2fall(riseheight >= thr_rise, 1));
                t_fall = t_bite(fall2rise(fallheight >= thr_fall, 1));
                
                [t_aligned, zscore_lift, zscore_lift_random] = AlignSignal2Event(fpdata_t, fpdata_zsignal, t_rise, pre, post);
                if ~isempty(Bite_events)
                    nrise = numel(t_rise);
                    firstbite = zeros(nrise, 1);
                    bitealign2rise = repmat(bite_timestamps, nrise, 1)-t_rise;
                    for k = 1:nrise
                        temp = bitealign2rise(k, :);
                        firstbite(k) = temp(find(temp >= 0, 1));
                    end
                    firstbite_all = [firstbite_all; firstbite];
                    [~, zscore_firstbite, zscore_firstbite_random] = AlignSignal2Event(fpdata_t, fpdata_zsignal, t_rise+firstbite, pre, post);
                end
                
                [~, zscore_fall, zscore_fall_random] = AlignSignal2Event(fpdata_t, fpdata_zsignal, t_fall, pre, post);
                if ~isempty(Bite_events)
                    nfall = numel(t_fall);
                    lastbite = zeros(nfall, 1);
                    bitealign2fall = repmat(bite_timestamps, nfall, 1)-t_fall;
                    for k = 1:nfall
                        temp = bitealign2fall(k, :);
                        lastbite(k) = temp(find(temp >= 0, 1));
                    end
                    lastbite_all = [lastbite_all; lastbite];
                    t_lastbite = t_fall+lastbite;
                    t_lastbite(lastbite > t_limit) = [];
                    [~, zscore_lastbite, zscore_lastbite_random] = AlignSignal2Event(fpdata_t, fpdata_zsignal, t_lastbite, pre, post);
                end
                
                zscore_lift_all = [zscore_lift_all zscore_lift];
                zscore_lift_random_all = [zscore_lift_random_all zscore_lift_random];
                zscore_fall_all = [zscore_fall_all zscore_fall];
                zscore_fall_random_all = [zscore_fall_random_all zscore_fall_random];
                if ~isempty(Bite_events)
                    zscore_firstbite_all = [zscore_firstbite_all zscore_firstbite];
                    zscore_firstbite_random_all = [zscore_firstbite_random_all zscore_firstbite_random];
                    zscore_lastbite_all = [zscore_lastbite_all zscore_lastbite];
                    zscore_lastbite_random_all = [zscore_lastbite_random_all zscore_lastbite_random];
                end
            end
        end
    end
end

t_cc = (shifttime:-1:-shifttime)'*ccstep/FrameRate;
for i = 1:nchannel
    figure;
    if ~isempty(Bite_events)
        subplot(2, 2, 1);
        yyaxis left;
        plot_gcamp_align2event(t_aligned, zscore_lift_random_all(:, :, i), zscore_lift_all(:, :, i), ['aligned to lift: ' lgdtext{i}]);
        set(gca, 'YColor', [0 0.5 0]);
        
        yyaxis right;
        h1 = cdfplot(firstbite_all);
        set(h1, 'Color', [1 0 0], 'LineWidth', 1.5);
        set(gca, 'GridLineStyle', ':', 'YColor', [1 0 0]);
        xlabel('Time (s)');
        ylabel('Cumulative Probability');
        title(['aligned to lift: ' lgdtext{i}]);
        
        subplot(2, 2, 2);
        plot_gcamp_align2event(t_aligned, zscore_firstbite_random_all(:, :, i), zscore_firstbite_all(:, :, i), ['aligned to first bite after lift: ' lgdtext{i}]);
        set(gca, 'ButtonDownFcn', @extract_figure);
        
        subplot(2, 2, 3);
        yyaxis left;
        plot_gcamp_align2event(t_aligned, zscore_fall_random_all(:, :, i), zscore_fall_all(:, :, i), ['aligned to fall: ' lgdtext{i}]);
        set(gca, 'YColor', [0 0.5 0]);
        
        yyaxis right;
        h1 = cdfplot(lastbite_all);
        set(h1, 'Color', [1 0 0], 'LineWidth', 1.5);
        line([t_limit t_limit], [0 1], 'Color', [1 0 0], 'LineStyle', '--', 'LineWidth', 1)
        set(gca, 'GridLineStyle', ':', 'YColor', [1 0 0]);
        xlabel('Time (s)');
        ylabel('Cumulative Probability');
        title(['aligned to fall: ' lgdtext{i}]);
        
        subplot(2, 2, 4);
        plot_gcamp_align2event(t_aligned, zscore_lastbite_random_all(:, :, i), zscore_lastbite_all(:, :, i), ['aligned to first bite after fall within ' num2str(t_limit) ' s : ' lgdtext{i}]);
        set(gca, 'ButtonDownFcn', @extract_figure);
    else
        subplot(1, 2, 1);
        plot_gcamp_align2event(t_aligned, zscore_lift_random_all(:, :, i), zscore_lift_all(:, :, i), ['aligned to lift: ' lgdtext{i}]);
        set(gca, 'ButtonDownFcn', @extract_figure);
        
        subplot(1, 2, 2);
        plot_gcamp_align2event(t_aligned, zscore_fall_random_all(:, :, i), zscore_fall_all(:, :, i), ['aligned to fall: ' lgdtext{i}]);
        set(gca, 'ButtonDownFcn', @extract_figure);
    end
    
    figure;
    subplot(1, 2, 1);
    plot_tj_individuals(t_cc, cc_all(:, :, i), [0.75 0.75 0.75], [0 0 0], 'Lag (s)', 'Correlation Coefficient', ['Z score with Z position: ' lgdtext{i}]);
    set(gca, 'ButtonDownFcn', @extract_figure);
    subplot(1, 2, 2);
    plot_tj_individuals(t_cc, cc2nose_all(:, :, i), [0.75 0.75 0.75], [0 0 0], 'Lag (s)', 'Correlation Coefficient', ['Z score with relative Z distance to nose: ' lgdtext{i}]);
    set(gca, 'ButtonDownFcn', @extract_figure);
    
    figure;
    hp(1) = plot_tj_MeanSEM(t_cc, cc_all(:, :, i), [0 1 1], [0 0 0.8], [], [], []);
    hp(2) = plot_tj_MeanSEM(t_cc, cc2nose_all(:, :, i), [1 0.41 0.16], [0.8 0 0], 'Lag (s)', 'Correlation Coefficient', lgdtext{i});
    yl = ylim;
    line([0 0], yl, 'Color', [0 0 0], 'LineStyle', '--', 'LineWidth', 1);
    legend(hp, {'Z score with Z position', 'Z score with Z distance to nose'});
    legend('boxoff');
end