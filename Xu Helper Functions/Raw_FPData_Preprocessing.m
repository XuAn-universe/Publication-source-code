function Raw_FPData_Preprocessing(signal_path, trigger_path, save_path, interval_exc, airPLS_check, lambda, globalbaseline, tbaseline, correct1st, lambda4trial1, smooth_check, swindow)
%*---------------------------------------------------------------------*
% Huang Lab
% Cold Spring Harbor Laboratory
% Author : Xu An, July 2020
% xan@cshl.edu
% Version: 1.0
%*---------------------------------------------------------------------*
% example:
% Raw_FPData_Preprocessing('D:\Fiber Photometry Data\050120_pellet_Fezf2CreER_Gcamp1\exp001\signal0.csv',...
%                          'D:\Fiber Photometry Data\050120_pellet_Fezf2CreER_Gcamp1\exp001\trigger0.csv',...
%                          'D:\Free-moving Feeding Data\050120_pellet_Fezf2CreER_Gcamp1', [680 720; 955 1200; 2422 inf], 1, 1e12, 0, 10.6, 1, 1e8);
if nargin < 12
    smooth_check = 0;
    swindow = 9;
    if nargin < 10
        correct1st = 0; % correct 1st trial from error in baseline correction
        lambda4trial1 = [];
        if nargin < 8
            globalbaseline = 1;
            tbaseline = []; % time period before door open for baseline estimation
            if nargin < 6
                lambda = 1e12;
                if nargin < 5
                    airPLS_check = 1;
                    if nargin < 4
                        interval_exc = [];
                    end
                end
            end
        end
    end
end
pre = 4; % for global baseline, keep 4 sec before trial starts; otherwise keep tbaseline sec before trial starts

parameter.interval_exc = interval_exc;
parameter.airPLS_check = airPLS_check;
parameter.lambda = lambda;
parameter.globalbaseline = globalbaseline;
parameter.tbaseline = tbaseline;
parameter.pre = pre;
parameter.correct1st = correct1st;
parameter.lambda4trial1 = lambda4trial1;
parameter.smooth_check = smooth_check;
parameter.swindow = swindow;

data = csvread(signal_path);
nchannel = size(data, 2) - 2;
if max(diff(data(:, 1))) > 1
    warndlg('There is frame/frames missing!', 'Warning');
    return
end
trigger = csvread(trigger_path);
ntrials = numel(trigger);
t = (data(:, 2)-trigger(1))/1000;
trigger = (trigger-trigger(1))/1000;

gcamp = data(1:2:end, 3:end);
isosbestic = data(2:2:end, 3:end);
tgcamp = t(1:2:end);
tisosbestic = t(2:2:end);
if globalbaseline
    index_exc = find(tgcamp+pre < 0);
else
    index_exc = find(tgcamp+tbaseline < 0);
end
gcamp(index_exc, :) = [];
isosbestic(index_exc, :) = [];
tgcamp(index_exc) = [];
tisosbestic(index_exc) = [];
if size(gcamp, 1) > size(isosbestic, 1)
    gcamp(end, :) = [];
    tgcamp(end, :) = [];
end
disp(['Sample Interval = ' num2str(mean(diff(tgcamp)))]);

gcamp_smoothed = NaN(size(gcamp));
isosbestic_smoothed = NaN(size(isosbestic));
gcamp_blcorrected = NaN(size(gcamp));
isosbestic_blcorrected = NaN(size(isosbestic));
gcamp_bl = NaN(size(gcamp));
isosbestic_bl = NaN(size(isosbestic));
zgcamp = NaN(size(gcamp));
zisosbestic = NaN(size(isosbestic));

zsignal = NaN(size(isosbestic));
zsignal_all = cell(nchannel, ntrials);

if ~isempty(interval_exc)
    if isfinite(interval_exc(end))
        index_inc = cell(1, size(interval_exc, 1)+1);
        for i = 1:numel(index_inc)
            if i ~= 1
                thead = interval_exc(i-1, 2);
            else
                thead = -inf;
            end
            if i ~= numel(index_inc)
                ttail = interval_exc(i, 1);
            else
                ttail = inf;
            end
            index_inc{i} = find(tgcamp < ttail & tgcamp > thead);
        end
    else
        index_inc = cell(1, size(interval_exc, 1));
        for i = 1:numel(index_inc)
            if i ~= 1
                thead = interval_exc(i-1, 2);
            else
                thead = -inf;
            end
            ttail = interval_exc(i, 1);
            index_inc{i} = find(tgcamp < ttail & tgcamp > thead);
        end
    end
end

for i = 1:nchannel
    % smoothing
    if smooth_check
        if isempty(interval_exc)
            gcamp_smoothed(:, i) = smooth(gcamp(:, i), swindow);
            isosbestic_smoothed(:, i) = smooth(isosbestic(:, i), swindow);
        else
            for j = 1:numel(index_inc)
                gcamp_smoothed(index_inc{j}, i) = smooth(gcamp(index_inc{j}, i), swindow);
                isosbestic_smoothed(index_inc{j}, i) = smooth(isosbestic(index_inc{j}, i), swindow);
            end
        end
        
        figure;
        subplot(2, 1, 1);
        plot(tgcamp, gcamp(:, i), '-k');
        hold on;
        plot(tgcamp, gcamp_smoothed(:, i), '-r');
        yl = ylim;
        line([trigger'; trigger'], repmat(yl', 1, ntrials), 'Color', 'g');
        xlim([tgcamp(1) tgcamp(end)]);
        xlabel('Time (s)');
        ylabel('F');
        title(['Gcamp (smooth): Channel ' num2str(i)]);
        set(gca, 'TickLength', [0 0], 'FontSize', 12, 'ButtonDownFcn', @extract_figure);
        
        subplot(2, 1, 2);
        plot(tisosbestic, isosbestic(:, i), '-k');
        hold on;
        plot(tisosbestic, isosbestic_smoothed(:, i), '-r');
        yl = ylim;
        line([trigger'; trigger'], repmat(yl', 1, ntrials), 'Color', 'g');
        xlim([tisosbestic(1) tisosbestic(end)]);
        xlabel('Time (s)');
        ylabel('F');
        title(['Isosbestic (smooth): Channel ' num2str(i)]);
        set(gca, 'TickLength', [0 0], 'FontSize', 12, 'ButtonDownFcn', @extract_figure);
        
        gcamp(:, i) = gcamp_smoothed(:, i);
        isosbestic(:, i) = isosbestic_smoothed(:, i);
    end
    
    % baseline correction
    if airPLS_check
        if isempty(interval_exc)
            [gcamp_blcorrected(:, i), gcamp_bl(:, i)] = airPLS(gcamp(:, i)', lambda);
            [isosbestic_blcorrected(:, i), isosbestic_bl(:, i)] = airPLS(isosbestic(:, i)', lambda);
        else
            for j = 1:numel(index_inc)
                [gcamp_blcorrected(index_inc{j}, i), gcamp_bl(index_inc{j}, i)] = airPLS(gcamp(index_inc{j}, i)', lambda);
                [isosbestic_blcorrected(index_inc{j}, i), isosbestic_bl(index_inc{j}, i)] = airPLS(isosbestic(index_inc{j}, i)', lambda);
            end
        end
        if correct1st
            % redo baseline correction for 1st trial
            if ~globalbaseline
                [gcamp_blcorrected(tgcamp < trigger(2)-tbaseline, i), gcamp_bl(tgcamp < trigger(2)-tbaseline, i)] = airPLS(gcamp(tgcamp < trigger(2)-tbaseline, i)', lambda4trial1);
                [isosbestic_blcorrected(tisosbestic < trigger(2)-tbaseline, i), isosbestic_bl(tisosbestic < trigger(2)-tbaseline, i)] = airPLS(isosbestic(tisosbestic < trigger(2)-tbaseline, i)', lambda4trial1);
            else
                [gcamp_blcorrected(tgcamp < trigger(2)-pre, i), gcamp_bl(tgcamp < trigger(2)-pre, i)] = airPLS(gcamp(tgcamp < trigger(2)-pre, i)', lambda4trial1);
                [isosbestic_blcorrected(tisosbestic < trigger(2)-pre, i), isosbestic_bl(tisosbestic < trigger(2)-pre, i)] = airPLS(isosbestic(tisosbestic < trigger(2)-pre, i)', lambda4trial1);
            end
        end
    else
        modelstr = 'y ~ b1*exp(-b2*x)+b3';
        opts = statset('fitnlm');
        opts.RobustWgtFun = []; % or 'bisquare'
        opts.MaxIter = 10000;
        if isempty(interval_exc)
            b0 = [0 0 0];
            mdl = fitnlm(tgcamp, gcamp(:, i), modelstr, b0, 'Options', opts);
            gcamp_bl(:, i) = predict(mdl, tgcamp);
            gcamp_blcorrected(:, i) = gcamp(:, i)-gcamp_bl(:, i);
            
            b0 = [0 0 0];
            mdl = fitnlm(tisosbestic, isosbestic(:, i), modelstr, b0, 'Options', opts);
            isosbestic_bl(:, i) = predict(mdl, tisosbestic);
            isosbestic_blcorrected(:, i) = isosbestic(:, i)-isosbestic_bl(:, i);
        else
            for j = 1:numel(index_inc)
                [gcamp_blcorrected(index_inc{j}, i), gcamp_bl(index_inc{j}, i)] = airPLS(gcamp(index_inc{j}, i)', lambda);
                [isosbestic_blcorrected(index_inc{j}, i), isosbestic_bl(index_inc{j}, i)] = airPLS(isosbestic(index_inc{j}, i)', lambda);
                
                b0 = [0 0 0];
                mdl = fitnlm(tgcamp(index_inc{j}), gcamp(index_inc{j}, i), modelstr, b0, 'Options', opts);
                gcamp_bl(index_inc{j}, i) = predict(mdl, tgcamp(index_inc{j}));
                gcamp_blcorrected(index_inc{j}, i) = gcamp(index_inc{j}, i)-gcamp_bl(index_inc{j}, i);
                
                b0 = [0 0 0];
                mdl = fitnlm(tisosbestic(index_inc{j}), isosbestic(index_inc{j}, i), modelstr, b0, 'Options', opts);
                isosbestic_bl(index_inc{j}, i) = predict(mdl, tisosbestic(index_inc{j}));
                isosbestic_blcorrected(index_inc{j}, i) = isosbestic(index_inc{j}, i)-isosbestic_bl(index_inc{j}, i);
            end
        end
        if correct1st
            % redo baseline correction for 1st trial
            if ~globalbaseline
                b0 = [0 0 0];
                mdl = fitnlm(tgcamp(tgcamp < trigger(2)-tbaseline), gcamp(tgcamp < trigger(2)-tbaseline, i), modelstr, b0, 'Options', opts);
                gcamp_bl(tgcamp < trigger(2)-tbaseline, i) = predict(mdl, tgcamp(tgcamp < trigger(2)-tbaseline));
                gcamp_blcorrected(tgcamp < trigger(2)-tbaseline, i) = gcamp(tgcamp < trigger(2)-tbaseline, i)-gcamp_bl(tgcamp < trigger(2)-tbaseline, i);
                
                b0 = [0 0 0];
                mdl = fitnlm(tisosbestic(tisosbestic < trigger(2)-tbaseline), isosbestic(tisosbestic < trigger(2)-tbaseline, i), modelstr, b0, 'Options', opts);
                isosbestic_bl(tisosbestic < trigger(2)-tbaseline, i) = predict(mdl, tisosbestic(tisosbestic < trigger(2)-tbaseline));
                isosbestic_blcorrected(tisosbestic < trigger(2)-tbaseline, i) = isosbestic(tisosbestic < trigger(2)-tbaseline, i)-isosbestic_bl(tisosbestic < trigger(2)-tbaseline, i);
            else
                b0 = [0 0 0];
                mdl = fitnlm(tgcamp(tgcamp < trigger(2)-pre), gcamp(tgcamp < trigger(2)-pre, i), modelstr, b0, 'Options', opts);
                gcamp_bl(tgcamp < trigger(2)-pre, i) = predict(mdl, tgcamp(tgcamp < trigger(2)-pre));
                gcamp_blcorrected(tgcamp < trigger(2)-pre, i) = gcamp(tgcamp < trigger(2)-pre, i)-gcamp_bl(tgcamp < trigger(2)-pre, i);
                
                b0 = [0 0 0];
                mdl = fitnlm(tisosbestic(tisosbestic < trigger(2)-pre), isosbestic(tisosbestic < trigger(2)-pre, i), modelstr, b0, 'Options', opts);
                isosbestic_bl(tisosbestic < trigger(2)-pre, i) = predict(mdl, tisosbestic(tisosbestic < trigger(2)-pre));
                isosbestic_blcorrected(tisosbestic < trigger(2)-pre, i) = isosbestic(tisosbestic < trigger(2)-pre, i)-isosbestic_bl(tisosbestic < trigger(2)-pre, i);
            end
        end
    end
    
    figure;
    subplot(2, 1, 1);
    plot(tgcamp, gcamp(:, i), '-k');
    hold on;
    plot(tgcamp, gcamp_bl(:, i), '-r');
    yl = ylim;
    line([trigger'; trigger'], repmat(yl', 1, ntrials), 'Color', 'g');
    xlim([tgcamp(1) tgcamp(end)]);
    xlabel('Time (s)');
    ylabel('F');
    title(['Gcamp (baseline correction): Channel ' num2str(i)]);
    set(gca, 'TickLength', [0 0], 'FontSize', 12, 'ButtonDownFcn', @extract_figure);
    
    subplot(2, 1, 2);
    plot(tisosbestic, isosbestic(:, i), '-k');
    hold on;
    plot(tisosbestic, isosbestic_bl(:, i), '-r');
    yl = ylim;
    line([trigger'; trigger'], repmat(yl', 1, ntrials), 'Color', 'g');
    xlim([tisosbestic(1) tisosbestic(end)]);
    xlabel('Time (s)');
    ylabel('F');
    title(['Isosbestic (baseline correction): Channel ' num2str(i)]);
    set(gca, 'TickLength', [0 0], 'FontSize', 12, 'ButtonDownFcn', @extract_figure);
    
    % z score
    if globalbaseline
        zgcamp(:, i) = (gcamp_blcorrected(:, i)-median(gcamp_blcorrected(:, i), 'omitnan'))/std(gcamp_blcorrected(:, i), 'omitnan');
        zisosbestic(:, i) = (isosbestic_blcorrected(:, i)-median(isosbestic_blcorrected(:, i), 'omitnan'))/std(isosbestic_blcorrected(:, i), 'omitnan');
    else
        zgcamp(:, i) = gcamp_blcorrected(:, i);
        zisosbestic(:, i) = isosbestic_blcorrected(:, i);
        for j = 1:ntrials
            if j ~= ntrials
                t_nexttrial = trigger(j+1);
            else
                t_nexttrial = inf;
            end
            gcamp_baseline = gcamp_blcorrected(tgcamp >= trigger(j)-tbaseline & tgcamp < trigger(j), i);
            gcamp_singletrial = gcamp_blcorrected(tgcamp >= trigger(j)-tbaseline & tgcamp < t_nexttrial-tbaseline, i);
            zgcamp_singletrial = (gcamp_singletrial-median(gcamp_baseline, 'omitnan'))/std(gcamp_baseline, 'omitnan');
            zgcamp(tgcamp >= trigger(j)-tbaseline & tgcamp < t_nexttrial-tbaseline, i) = zgcamp_singletrial;
            isosbestic_baseline = isosbestic_blcorrected(tisosbestic >= trigger(j)-tbaseline & tisosbestic < trigger(j), i);
            isosbestic_singletrial = isosbestic_blcorrected(tisosbestic >= trigger(j)-tbaseline & tisosbestic < t_nexttrial-tbaseline, i);
            zisosbestic_singletrial = (isosbestic_singletrial-median(isosbestic_baseline, 'omitnan'))/std(isosbestic_baseline, 'omitnan');
            zisosbestic(tisosbestic >= trigger(j)-tbaseline & tisosbestic < t_nexttrial-tbaseline, i) = zisosbestic_singletrial;
        end
        
    end
    
    figure;
    subplot(4, 1, 1);
    plot(tgcamp, zgcamp(:, i), '-k');
    hold on;
    yl = ylim;
    line([trigger'; trigger'], repmat(yl', 1, ntrials), 'Color', 'g');
    xlim([tgcamp(1) tgcamp(end)]);
    xlabel('Time (s)');
    ylabel('Z-score');
    title(['Gcamp (Z-score): Channel ' num2str(i)]);
    set(gca, 'TickLength', [0 0], 'FontSize', 12, 'ButtonDownFcn', @extract_figure);
    
    subplot(4, 1, 2);
    plot(tisosbestic, zisosbestic(:, i), '-k');
    hold on;
    yl = ylim;
    line([trigger'; trigger'], repmat(yl', 1, ntrials), 'Color', 'g');
    xlim([tisosbestic(1) tisosbestic(end)]);
    xlabel('Time (s)');
    ylabel('Z-score');
    title(['Isosbestic (Z-score): Channel ' num2str(i)]);
    set(gca, 'TickLength', [0 0], 'FontSize', 12, 'ButtonDownFcn', @extract_figure);
    
    % isosbestic correction
    index_temp = find(~isnan(zgcamp(:, i)));
    mdl = fitlm(zisosbestic(index_temp, i), zgcamp(index_temp, i), 'RobustOpts', 'on');
    fitzisosbestic = predict(mdl, zisosbestic(:, i));
    zsignal(:, i) = zgcamp(:, i) - fitzisosbestic;
    
    subplot(4, 1, 3);
    plot(tgcamp, zgcamp(:, i), '-k');
    hold on;
    plot(tgcamp, fitzisosbestic, '-r');
    yl = ylim;
    line([trigger'; trigger'], repmat(yl', 1, ntrials), 'Color', 'g');
    xlim([tgcamp(1) tgcamp(end)]);
    xlabel('Time (s)');
    ylabel('Z-score');
    title(['Isosbestic fit to Gcamp: Channel ' num2str(i)]);
    set(gca, 'TickLength', [0 0], 'FontSize', 12, 'ButtonDownFcn', @extract_figure);
    
    subplot(4, 1, 4);
    plot(tgcamp, zsignal(:, i), '-k');
    hold on;
    line([tgcamp(1) tgcamp(end)], [0 0], 'Color', 'r');
    yl = ylim;
    line([trigger'; trigger'], repmat(yl', 1, ntrials), 'Color', 'g');
    xlim([tisosbestic(1) tisosbestic(end)]);
    xlabel('Time (s)');
    ylabel('Z-score');
    title(['Signal: Channel ' num2str(i)]);
    set(gca, 'TickLength', [0 0], 'FontSize', 12, 'ButtonDownFcn', @extract_figure);
    
    for j = 1:ntrials
        if j ~= ntrials
            t_nexttrial = trigger(j+1);
        else
            t_nexttrial = inf;
        end
        if globalbaseline
            t_trial = tgcamp(tgcamp >= trigger(j)-pre & tgcamp < t_nexttrial-pre)-trigger(j);
            zsignal_trial = zsignal(tgcamp >= trigger(j)-pre & tgcamp < t_nexttrial-pre, i);
            zsignal_all{i, j} = [t_trial zsignal_trial];
        else
            t_trial = tgcamp(tgcamp >= trigger(j)-tbaseline & tgcamp < t_nexttrial-tbaseline)-trigger(j);
            zsignal_trial = zsignal(tgcamp >= trigger(j)-tbaseline & tgcamp < t_nexttrial-tbaseline, i);
            zsignal_all{i, j} = [t_trial zsignal_trial];
        end
    end
end
save([save_path '\FPData.mat'], 'zsignal_all', 'parameter');