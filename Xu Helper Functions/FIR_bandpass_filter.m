function filter_signal = FIR_bandpass_filter(signal, SampleRate, passband, transition_width, nan_window, check_filter)
% FIR filter design
nyquist = SampleRate/2; % theoretical limit
pband1 = passband(1);
pband2 = passband(2);
filterlen = round(1/pband1*SampleRate*1); % larger the order, better the frequency response
if numel(signal) <= 3*filterlen
    filterlen = round(2*SampleRate);
end

% transition_width = 0.1; % usually 0.1~0.25, the sharpeness of the edge

ffrequencies   = [0 (1-transition_width)*pband1 pband1 pband2 (1+transition_width)*pband2 nyquist]/nyquist;
idealresponse  = [0 0 1 1 0 0]; % Band-pass filter
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

filter_signal = filtfilt(filterweights, 1, fillmissing(signal, 'movmean', nan_window)); % apply filter to the data