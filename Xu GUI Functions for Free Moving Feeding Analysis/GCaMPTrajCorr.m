function GCaMPTrajCorr(app, Exp_Path, FrameRate, nmedian)
%*---------------------------------------------------------------------*
% Huang Lab
% Duke University
% Author : Xu An, Nov 2022
% xu.an@duke.edu
% Version: 1.0
%*---------------------------------------------------------------------*
value = app.TrialsListBox.Value;
value = sort(value);
trials = numel(value);

if app.TrackingDataCheckBox.Value
    label_table = table2array(app.UITable.Data);
else
    helpdlg('''Tracking Data'' needs to be checked');
    return;
end

check_filter = 0;

ccstep = 1; % number of data point, default 1, shift 1 data point each time
shifttime = 120; % half of shift time, default 120, shift 1 sec for both left and right sides
cc_all = [];
trialID = cell(0);

load([Exp_Path '\Analysis_Session.mat'], 'Video_annotation');
    
Bite_events = [];
% get bite events
try
    audiolocation = Exp_Path(1:end-7);
    temp = load([audiolocation '\Detected_Bite_Events.mat']);
    Bite_events = temp.Audio_analysis;
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

for i = 1:trials
    if i == 1
        % FIR filter design
        nyquist = FrameRate/2; % theoretical limit
        pband = 5;
        filterlen = round(15/pband*FrameRate*1); % larger the order, better the frequency response
        
        transition_width = 0.2; % usually 0.1~0.25, the sharpeness of the edge
        
        ffrequencies   = [0 pband (1+transition_width)*pband nyquist]/nyquist;
        idealresponse  = [0 0 1 1]; % High-pass filter
        filterweights = firls(filterlen, ffrequencies, idealresponse);
        
        % Check the Filter Quality
        if check_filter
            filterweightsW = zscore(filterweights); % with z-score
            plot(ffrequencies*nyquist, idealresponse, 'r');
            hold on;
            
            fft_filtkern  = abs(fft(filterweightsW));
            fft_filtkern  = fft_filtkern./max(fft_filtkern); % normalized to 1.0 for visual comparison ease
            
            hz_filtkern = linspace(0, nyquist, floor(length(fft_filtkern)/2+1));
            plot(hz_filtkern, fft_filtkern(1:ceil(length(fft_filtkern)/2)), 'b');
            set(gca, 'ylim', [-.1 1.1], 'xlim', [0 nyquist]);
            set(gca, 'FontSize', 12, 'ButtonDownFcn', @extract_figure);
            legend({'Ideal', 'Best fit'});
            freqsidx = dsearchn(hz_filtkern', ffrequencies'*nyquist);
            title([ 'SSE: ' num2str(sum((idealresponse-fft_filtkern(freqsidx)).^2 )) ]);
        end
    end

    if ~Video_annotation(value(i)).Disgard
        bite_timestamps = [];
        if ~isempty(Bite_events)
            bite_timestamps = Bite_events(value(i)).time_bites;
            bite_timestamps = sort(bite_timestamps);
        else
            if ~isempty(Video_annotation(value(i)).time_ready2bite) && ~isempty(Video_annotation(value(i)).time_feeding_end)
                bite_timestamps(1) = Video_annotation(value(i)).time_ready2bite;
                bite_timestamps(2) = Video_annotation(value(i)).time_feeding_end;
            end
        end
        
        fpdata = fpdata_all.zsignal_all(:, value(i));
        fpdata_zsignal = [];
        fpdata_t = [];
        for j = 1:nchannel
            fpdata_zsignal(:, j) = fpdata{j}(:, 2);
        end
        fpdata_t = fpdata{1}(:, 1);
        SampleRate = mean(diff(fpdata_t));
        
        [~, ~, ~, ~, ~, ~, ~, znose, ~, ~, ~, ~] = trajectory_postprocessing(3, Exp_Path, value(i), label_table, nmedian, FrameRate);
        [~, ~, ~, ~, ~, ~, ~, zpawl, ~, ~, ~, ~] = trajectory_postprocessing(10, Exp_Path, value(i), label_table, nmedian, FrameRate);
        [~, ~, ~, ~, ~, ~, ~, zpawr, ~, ~, ~, ~] = trajectory_postprocessing(16, Exp_Path, value(i), label_table, nmedian, FrameRate);
        paw2nose = (zpawl+zpawr)/2-znose;
        paw2nose = paw2nose(:, 1);
        t = (1:size(paw2nose, 1))'/FrameRate;
        
        if numel(bite_timestamps) >= 2
            paw2nose_ROI = paw2nose(t >= bite_timestamps(1)-shifttime/FrameRate & t <= bite_timestamps(end)+shifttime/FrameRate);
            paw2nose_ROI = fillmissing(paw2nose_ROI, 'spline');
%             paw2nose_ROI = fillmissing(paw2nose_ROI, 'movmean', 12);
            paw2nose_ROI = paw2nose_ROI-filtfilt(filterweights, 1, paw2nose_ROI); % apply filter to the data
            paw2nose(t >= bite_timestamps(1)-shifttime/FrameRate & t <= bite_timestamps(end)+shifttime/FrameRate) = paw2nose_ROI;
            
            t_bite = t(t >= bite_timestamps(1) & t <= bite_timestamps(end));
            fpdata_zsignal_bite = fpdata_zsignal(fpdata_t >= bite_timestamps(1) & fpdata_t <= bite_timestamps(end), :);
            fpdata_t_bite = fpdata_t(fpdata_t >= bite_timestamps(1) & fpdata_t <= bite_timestamps(end));
            temp = repmat(t_bite', numel(fpdata_t_bite), 1)-fpdata_t_bite;
            [~, index_match] = min(abs(temp), [], 2);
            % compute correlation coefficients
            cc = zeros(2*shifttime+1, 1, nchannel);
            for j = -shifttime:1:shifttime
                paw2nose_shift = circshift(paw2nose, j*ccstep);
                paw2nose_bite_shift = paw2nose_shift(t >= bite_timestamps(1) & t <= bite_timestamps(end));
                paw2nose_bite_match = paw2nose_bite_shift(index_match);
                for k = 1:nchannel
                    R = corrcoef(paw2nose_bite_match, fpdata_zsignal_bite(:, k), 'Rows', 'pairwise');
                    cc(j+shifttime+1, 1, k) = R(2);
                end
                
                if j == 0
                    figure;
                    yyaxis left;
                    plot(t_bite, paw2nose_bite_shift, '-k');
                    ylabel('Hand-to-nose distance');
                    yyaxis right;
                    plot(fpdata_t_bite, fpdata_zsignal_bite(:, 2), '-g');
                    set(gca, 'TickDir', 'out', 'FontSize', 12);
                    box off;
                    xlabel('Time (s)');
                    ylabel('Z-score (Channel 2)');
                    title(['Trial ' num2str(value(i))]);
                end
            end
            cc_all = [cc_all cc];
            trialID(end+1) = {num2str(value(i))};
        end
    end
end

t_cc = (-shifttime:1:shifttime)'*ccstep/FrameRate;
for i = 1:nchannel
    figure;
    plot_tj_individuals(t_cc, cc_all(:, :, i), [0.75 0.75 0.75], [0 0 0], 'Lag (s)', 'Correlation coefficient', ['Z-score with hand-to-nose distance: ' lgdtext{i}], trialID);
    
    figure;
    hp = plot_tj_MeanSEM(t_cc, cc_all(:, :, i), [1 0.41 0.16], [0.8 0 0], 'Lag (s)', 'Correlation coefficient', lgdtext{i});
    yl = ylim;
    line([0 0], yl, 'Color', [0 0 0], 'LineStyle', '--', 'LineWidth', 1);
    legend(hp, 'Z-score with hand-to-nose distance');
    legend('boxoff');
    
    cc_avg = mean(cc_all(:, :, i), 2);
    [cc_peak, maxID] = max(cc_avg);
    t_peak = t_cc(maxID);
    disp([lgdtext{i} ': peak cc = ' num2str(cc_peak) '; peak time = ' num2str(t_peak)]);
end

result.cc_all = cc_all;
assignin('base', 'result', result);