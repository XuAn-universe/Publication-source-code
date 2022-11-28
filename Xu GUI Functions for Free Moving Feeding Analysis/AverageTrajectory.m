function AverageTrajectory(app, Exp_Path, FrameRate, nmedian)
persistent path
filtering4phase = 0; % whether to filter z trajectory to compute phase
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

n_NI = 0;
x_all_NI = cell(0);
y_all_NI = cell(0);
z_all_NI = cell(0);
speed_all_NI = cell(0);
acceleration_all_NI = cell(0);
frames_NI = 0;

bite_firstlast = [];

n_I = 0;
x_all_I = cell(0);
y_all_I = cell(0);
z_all_I = cell(0);
speed_all_I = cell(0);
acceleration_all_I = cell(0);
frames_I = 0;
laserstart_all = [];
laserstop_all = [];

zbite_NI = [];
zspeedbite_NI = [];
xyzspeedbite_NI = [];
zaccelerationbite_NI = [];
xyzaccelerationbite_NI = [];
zphasebite_NI =[];

zbite_I = [];
zspeedbite_I = [];
xyzspeedbite_I = [];
zaccelerationbite_I = [];
xyzaccelerationbite_I = [];
zphasebite_I =[];

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
trange = app.bitesEditField.Value; % time range for bite alignment
trange = trange*FrameRate;

nbite_NI = zeros(1, ninhibition);
nbite_I = zeros(1, ninhibition);

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
    
    hz_filtkern = linspace(0, nyquist, floor(length(fft_filtkern)/2+1));
    plot(hz_filtkern,fft_filtkern(1:ceil(length(fft_filtkern)/2)),'b');
    
    set(gca,'ylim',[-.1 1.1],'xlim',[0 nyquist]);
    set(gca, 'TickLength', [0 0], 'FontSize', 12);
    legend({'ideal';'best fit'});
    
    freqsidx = dsearchn(hz_filtkern',ffrequencies'*nyquist);
    title([ 'SSE: ' num2str(sum((idealresponse-fft_filtkern(freqsidx)).^2 )) ]);
end

for i = 1:numel(Exp_Paths)
    load([Exp_Paths{i} '\Analysis_Session.mat'], 'Video_annotation');
    try
        audiolocation = Exp_Paths{i}(1:end-7);
        temp = load([audiolocation '\Detected_Bite_Events.mat']);
        Bite_events = temp.Audio_analysis;
    catch
        errordlg('Please detect the bites first.', 'Error');
        return;
    end
    for j = 1:numel(Video_annotation)
        if ~Video_annotation(j).Disgard
            events_related2drop = [];
            try
                temp = load([Exp_Paths{i} '\LabelledEvents' num2str(j) '.mat']);
                LabelledEvents = temp.LabelledEvents;
                MouthOpen = LabelledEvents.MouthOpen;
                TongueOut = LabelledEvents.TongueOut;
                PawLReachStart = LabelledEvents.PawLReachStart;
                PawRReachStart = LabelledEvents.PawRReachStart;
                events_related2drop = [MouthOpen; TongueOut; PawLReachStart; PawRReachStart];
            catch
                errordlg([Exp_Paths{i} '\LabelledEvents' num2str(j) '.mat is missing!'], 'Error');
            end
            
            if app.sDelay4sOnButton.Value
                if any(events_related2drop >= 4 & events_related2drop <= 8) % mouse drops pasta during inhibition
                    continue;
                end
            end
            
            [~, ~, ~, ~, ~, x, y, z, speed, acceleration, laserstart, laserstop] = trajectory_postprocessing(pointID, Exp_Paths{i}, j,...
                label_table, nmedian, FrameRate);
            t = (1:size(x, 1))'/FrameRate;
            
            bite_timestamps = Bite_events(j).time_bites;
            if app.NoLightButton.Value
                laserstart = [];
                laserstop = [];
            else
                laser_timestamps = Bite_events(j).laser_timestamps;
            end
            
            % compute the phase below
            if app.sDelay4sOnButton.Value
                zphase = nan(size(z, 1), 1);
                if numel(bite_timestamps) >= 2
                    zshort = z(t >= bite_timestamps(1) & t <= bite_timestamps(end), 1);
                    try
                        if filtering4phase
                            zfiltered = filtfilt(filterweights, 1, fillmissing(zshort, 'movmedian', nmedian)); % apply bandpass filter to the data
                            zhilbert = hilbert(zfiltered);
                        else
                            zhilbert = hilbert(fillmissing(zshort, 'movmedian', nmedian));
                        end
                        zphase(t >= bite_timestamps(1) & t <= bite_timestamps(end)) = angle(zhilbert)/pi*180; % this angle is the cosine angle
                    end
                end
            end

            if ~isempty(laserstart)
                n_I = n_I+1;
                if numel(laserstart) ~= ninhibition
                    laserstart_all = [laserstart_all round((Bite_events(j).laser_timestamps(1, :))'*FrameRate)];
                    laserstop_all = [laserstop_all round((Bite_events(j).laser_timestamps(2, :))'*FrameRate)];
                else
                    laserstart_all = [laserstart_all laserstart];
                    laserstop_all = [laserstop_all laserstop];
                end
                
                x_all_I{n_I} = x(:, 1);
                y_all_I{n_I} = y(:, 1);
                z_all_I{n_I} = z(:, 1);
                speed_all_I{n_I} = speed;
                acceleration_all_I{n_I} = acceleration;
                frames_I = max(frames_I, size(x, 1));
                
                if ~isempty(bite_timestamps)
                    for k = 1:ninhibition
                        for m = 1:numel(bite_timestamps)
                            timestamp = bite_timestamps(m);
                            if timestamp >= laser_timestamps(2*k-1) && timestamp <= laser_timestamps(2*k)
                                [~, id] = min(abs(t-timestamp));
                                if id-trange >= 1
                                    nbite_I(k) = nbite_I(k)+1;
                                    zbite_I(:, nbite_I(k), k) = z(id-trange:id+trange, 1);
                                    zspeedbite_I(:, nbite_I(k), k) = abs(speed(id-trange:id+trange, 3));
                                    xyzspeedbite_I(:, nbite_I(k), k) = speed(id-trange:id+trange, 7);
                                    zaccelerationbite_I(:, nbite_I(k), k) = abs(acceleration(id-trange:id+trange, 3));
                                    xyzaccelerationbite_I(:, nbite_I(k), k) = abs(acceleration(id-trange:id+trange, 7));
                                    
                                    zphasebite_I(:, nbite_I(k), k) = zphase(id-trange:id+trange, 1);
                                end
                            end
                        end
                    end
                end
            else
                n_NI = n_NI+1;
                x_all_NI{n_NI} = x(:, 1);
                y_all_NI{n_NI} = y(:, 1);
                z_all_NI{n_NI} = z(:, 1);
                speed_all_NI{n_NI} = speed;
                acceleration_all_NI{n_NI} = acceleration;
                frames_NI = max(frames_NI, size(x, 1));
                
                if numel(bite_timestamps) >= 2
                    bite_firstlast(n_NI, :) = [bite_timestamps(1) bite_timestamps(end)];
                else
                    bite_firstlast(n_NI, :) = [NaN NaN];
                end
                
                if ~isempty(bite_timestamps)
                    for k = 1:ninhibition
                        for m = 1:numel(bite_timestamps)
                            timestamp = bite_timestamps(m);
                            if timestamp >= tinhibition(2*k-1) && timestamp <= tinhibition(2*k)
                                [~, id] = min(abs(t-timestamp));
                                if id-trange >= 1
                                    nbite_NI(k) = nbite_NI(k)+1;
                                    zbite_NI(:, nbite_NI(k), k) = z(id-trange:id+trange, 1);
                                    zspeedbite_NI(:, nbite_NI(k), k) = abs(speed(id-trange:id+trange, 3));
                                    xyzspeedbite_NI(:, nbite_NI(k), k) = speed(id-trange:id+trange, 7);
                                    zaccelerationbite_NI(:, nbite_NI(k), k) = abs(acceleration(id-trange:id+trange, 3));
                                    xyzaccelerationbite_NI(:, nbite_NI(k), k) = abs(acceleration(id-trange:id+trange, 7));
                                    
                                    if app.sDelay4sOnButton.Value
                                        zphasebite_NI(:, nbite_NI(k), k) = zphase(id-trange:id+trange, 1);
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end

x_matrix_NI = NaN(frames_NI, n_NI);
y_matrix_NI = NaN(frames_NI, n_NI);
z_matrix_NI = NaN(frames_NI, n_NI);
speed_matrix_NI = NaN(frames_NI, 7, n_NI);
acceleration_matrix_NI = NaN(frames_NI, 7, n_NI);

for i = 1:ninhibition
    zlight_all_NI{i} = [];
    zspeedlight_all_NI{i} = [];
    xyzspeedlight_all_NI{i} = [];
    zaccelerationlight_all_NI{i} = [];
    xyzaccelerationlight_all_NI{i} = [];
    
    zlight_mean_NI{i} = [];
    zspeedlight_mean_NI{i} = [];
    xyzspeedlight_mean_NI{i} = [];
    zaccelerationlight_mean_NI{i} = [];
    xyzaccelerationlight_mean_NI{i} = [];
    
    zlight_var_NI{i} = [];
    zspeedlight_var_NI{i} = [];
    xyzspeedlight_var_NI{i} = [];
    zaccelerationlight_var_NI{i} = [];
    xyzaccelerationlight_var_NI{i} = [];
end
for i = 1:n_NI
    frames_temp = numel(x_all_NI{i});
    x_matrix_NI(1:frames_temp, i) = x_all_NI{i};
    y_matrix_NI(1:frames_temp, i) = y_all_NI{i};
    z_matrix_NI(1:frames_temp, i) = z_all_NI{i};
    speed_matrix_NI(1:frames_temp, :, i) = abs(speed_all_NI{i});
    acceleration_matrix_NI(1:frames_temp, :, i) = abs(acceleration_all_NI{i});
    
    t = (1:frames_temp)'/FrameRate;
    for j = 1:ninhibition
        zlight_all_NI{j} = [zlight_all_NI{j}; z_all_NI{i}(t >= tinhibition(2*j-1) & t <= tinhibition(2*j))];
        zspeedlight_all_NI{j} = [zspeedlight_all_NI{j}; abs(speed_all_NI{i}(t >= tinhibition(2*j-1) & t <= tinhibition(2*j), 3))];
        xyzspeedlight_all_NI{j} = [xyzspeedlight_all_NI{j}; speed_all_NI{i}(t >= tinhibition(2*j-1) & t <= tinhibition(2*j), 7)];
        zaccelerationlight_all_NI{j} = [zaccelerationlight_all_NI{j}; abs(acceleration_all_NI{i}(t >= tinhibition(2*j-1) & t <= tinhibition(2*j), 3))];
        xyzaccelerationlight_all_NI{j} = [xyzaccelerationlight_all_NI{j}; abs(acceleration_all_NI{i}(t >= tinhibition(2*j-1) & t <= tinhibition(2*j), 7))];
        
        zlight_mean_NI{j} = [zlight_mean_NI{j}; mean(z_all_NI{i}(t >= tinhibition(2*j-1) & t <= tinhibition(2*j)), 'omitnan')];
        zspeedlight_mean_NI{j} = [zspeedlight_mean_NI{j}; mean(abs(speed_all_NI{i}(t >= tinhibition(2*j-1) & t <= tinhibition(2*j), 3)), 'omitnan')];
        xyzspeedlight_mean_NI{j} = [xyzspeedlight_mean_NI{j}; mean(speed_all_NI{i}(t >= tinhibition(2*j-1) & t <= tinhibition(2*j), 7), 'omitnan')];
        zaccelerationlight_mean_NI{j} = [zaccelerationlight_mean_NI{j}; mean(abs(acceleration_all_NI{i}(t >= tinhibition(2*j-1) & t <= tinhibition(2*j), 3)), 'omitnan')];
        xyzaccelerationlight_mean_NI{j} = [xyzaccelerationlight_mean_NI{j}; mean(abs(acceleration_all_NI{i}(t >= tinhibition(2*j-1) & t <= tinhibition(2*j), 7)), 'omitnan')];
        
        zlight_var_NI{j} = [zlight_var_NI{j}; var(z_all_NI{i}(t >= tinhibition(2*j-1) & t <= tinhibition(2*j)), 'omitnan')];
        zspeedlight_var_NI{j} = [zspeedlight_var_NI{j}; var(abs(speed_all_NI{i}(t >= tinhibition(2*j-1) & t <= tinhibition(2*j), 3)), 'omitnan')];
        xyzspeedlight_var_NI{j} = [xyzspeedlight_var_NI{j}; var(speed_all_NI{i}(t >= tinhibition(2*j-1) & t <= tinhibition(2*j), 7), 'omitnan')];
        zaccelerationlight_var_NI{j} = [zaccelerationlight_var_NI{j}; var(abs(acceleration_all_NI{i}(t >= tinhibition(2*j-1) & t <= tinhibition(2*j), 3)), 'omitnan')];
        xyzaccelerationlight_var_NI{j} = [xyzaccelerationlight_var_NI{j}; var(abs(acceleration_all_NI{i}(t >= tinhibition(2*j-1) & t <= tinhibition(2*j), 7)), 'omitnan')];
    end
end
t_NI = (1:frames_NI)'/FrameRate;
time_range = [app.timerangesEditField.Value app.toEditField.Value];

if ~app.NoLightButton.Value
    x_matrix_I = NaN(frames_I, n_I);
    y_matrix_I = NaN(frames_I, n_I);
    z_matrix_I = NaN(frames_I, n_I);
    speed_matrix_I = NaN(frames_I, 7, n_I);
    acceleration_matrix_I = NaN(frames_I, 7, n_I);
    
    for i = 1:ninhibition
        zlight_all_I{i} = [];
        zspeedlight_all_I{i} = [];
        xyzspeedlight_all_I{i} = [];
        zaccelerationlight_all_I{i} = [];
        xyzaccelerationlight_all_I{i} = [];
        
        zlight_mean_I{i} = [];
        zspeedlight_mean_I{i} = [];
        xyzspeedlight_mean_I{i} = [];
        zaccelerationlight_mean_I{i} = [];
        xyzaccelerationlight_mean_I{i} = [];
        
        zlight_var_I{i} = [];
        zspeedlight_var_I{i} = [];
        xyzspeedlight_var_I{i} = [];
        zaccelerationlight_var_I{i} = [];
        xyzaccelerationlight_var_I{i} = [];
    end
    for i = 1:n_I
        frames_temp = numel(x_all_I{i});
        x_matrix_I(1:frames_temp, i) = x_all_I{i};
        y_matrix_I(1:frames_temp, i) = y_all_I{i};
        z_matrix_I(1:frames_temp, i) = z_all_I{i};
        speed_matrix_I(1:frames_temp, :, i) = abs(speed_all_I{i});
        acceleration_matrix_I(1:frames_temp, :, i) = abs(acceleration_all_I{i});
        
        for j = 1:ninhibition
            zlight_all_I{j} = [zlight_all_I{j}; z_all_I{i}(laserstart_all(j, i):laserstop_all(j, i))];
            zspeedlight_all_I{j} = [zspeedlight_all_I{j}; abs(speed_all_I{i}(laserstart_all(j, i):laserstop_all(j, i), 3))];
            xyzspeedlight_all_I{j} = [xyzspeedlight_all_I{j}; speed_all_I{i}(laserstart_all(j, i):laserstop_all(j, i), 7)];
            zaccelerationlight_all_I{j} = [zaccelerationlight_all_I{j}; abs(acceleration_all_I{i}(laserstart_all(j, i):laserstop_all(j, i), 3))];
            xyzaccelerationlight_all_I{j} = [xyzaccelerationlight_all_I{j}; abs(acceleration_all_I{i}(laserstart_all(j, i):laserstop_all(j, i), 7))];
            
            zlight_mean_I{j} = [zlight_mean_I{j}; mean(z_all_I{i}(laserstart_all(j, i):laserstop_all(j, i)), 'omitnan')];
            zspeedlight_mean_I{j} = [zspeedlight_mean_I{j}; mean(abs(speed_all_I{i}(laserstart_all(j, i):laserstop_all(j, i), 3)), 'omitnan')];
            xyzspeedlight_mean_I{j} = [xyzspeedlight_mean_I{j}; mean(speed_all_I{i}(laserstart_all(j, i):laserstop_all(j, i), 7), 'omitnan')];
            zaccelerationlight_mean_I{j} = [zaccelerationlight_mean_I{j}; mean(abs(acceleration_all_I{i}(laserstart_all(j, i):laserstop_all(j, i), 3)), 'omitnan')];
            xyzaccelerationlight_mean_I{j} = [xyzaccelerationlight_mean_I{j}; mean(abs(acceleration_all_I{i}(laserstart_all(j, i):laserstop_all(j, i), 7)), 'omitnan')];
            
            zlight_var_I{j} = [zlight_var_I{j}; var(z_all_I{i}(laserstart_all(j, i):laserstop_all(j, i)), 'omitnan')];
            zspeedlight_var_I{j} = [zspeedlight_var_I{j}; var(abs(speed_all_I{i}(laserstart_all(j, i):laserstop_all(j, i), 3)), 'omitnan')];
            xyzspeedlight_var_I{j} = [xyzspeedlight_var_I{j}; var(speed_all_I{i}(laserstart_all(j, i):laserstop_all(j, i), 7), 'omitnan')];
            zaccelerationlight_var_I{j} = [zaccelerationlight_var_I{j}; var(abs(acceleration_all_I{i}(laserstart_all(j, i):laserstop_all(j, i), 3)), 'omitnan')];
            xyzaccelerationlight_var_I{j} = [xyzaccelerationlight_var_I{j}; var(abs(acceleration_all_I{i}(laserstart_all(j, i):laserstop_all(j, i), 7)), 'omitnan')];
        end
    end
    laserstart = mean(laserstart_all, 2)/FrameRate;
    laserstop = mean(laserstop_all, 2)/FrameRate;
    t_I = (1:frames_I)'/FrameRate;
    
%     plot_tj_multitrial(time_range, t_NI, x_matrix_NI, t_I, x_matrix_I, laserstart, laserstop, 'X (mm)', []);
%     plot_tj_multitrial(time_range, t_NI, y_matrix_NI, t_I, y_matrix_I, laserstart, laserstop, 'Y (mm)', []);
    figure;
    plot_tj_multitrial(time_range, t_NI, z_matrix_NI, t_I, z_matrix_I, laserstart, laserstop, 'Z (mm)', []);
    title_text = {'X', 'Y', 'Z', 'XY', 'XZ', 'YZ', 'XYZ'};
    for i = [3 7]
        figure;
        plot_tj_multitrial(time_range, t_NI, squeeze(speed_matrix_NI(:, i, :)), t_I, squeeze(speed_matrix_I(:, i, :)),...
            laserstart, laserstop, 'Speed (mm/s)', title_text{i});
        figure;
        plot_tj_multitrial(time_range, t_NI, squeeze(acceleration_matrix_NI(:, i, :)), t_I, squeeze(acceleration_matrix_I(:, i, :)),...
            laserstart, laserstop, 'Absolute Acceleration (mm/s^{2})', title_text{i});
    end
else
%     plot_tj_multitrial(time_range, t_NI, x_matrix_NI, [], [], [], [], 'X (mm)', []);
%     plot_tj_multitrial(time_range, t_NI, y_matrix_NI, [], [], [], [], 'Y (mm)', []);
    figure;
    plot_tj_multitrial(time_range, t_NI, z_matrix_NI, [], [], [], [], 'Z (mm)', []);
    title_text = {'X', 'Y', 'Z', 'XY', 'XZ', 'YZ', 'XYZ'};
    for i = [3 7]
        figure;
        plot_tj_multitrial(time_range, t_NI, squeeze(speed_matrix_NI(:, i, :)), [], [],...
            [], [], 'Speed (mm/s)', title_text{i});
        figure;
        plot_tj_multitrial(time_range, t_NI, squeeze(acceleration_matrix_NI(:, i, :)), [], [],...
            [], [], 'Absolute Acceleration (mm/s^{2})', title_text{i});
    end
end
try
    % Fourier transform of Z
    [~, temp] = max(diff(bite_firstlast, 1, 2));
    maxz = sum((t_NI >= bite_firstlast(temp, 1) & t_NI <= bite_firstlast(temp, 2)));
    n = 0;
    Fn = FrameRate/2;
    for i = 1:n_NI
        z = z_all_NI{i};
        t = (1:numel(z))'/FrameRate;
        z(t < bite_firstlast(i, 1) | t > bite_firstlast(i, 2)) = [];
        z = fillmissing(z, 'movmedian', nmedian);
        if ~any(isnan(z))
            n = n+1;
            winfunc = hanning(numel(z));
            NFFT = 2^nextpow2(maxz);         % Next highest power of 2 greater than length(x).
            NumUniquePts = ceil((NFFT+1)/2);
            f = (0:NumUniquePts-1)*2*Fn/NFFT;
            signal = (z-mean(z)).*winfunc;
            FFTX = fft(signal, NFFT);               % Take FFT, padding with zeros. length(FFTX)==NFFT
            FFTX = FFTX(1:NumUniquePts);            % FFT is symmetric, throw away second half
            MX = abs(FFTX);                         % Take magnitude of X, also equal to sqrt(FFTX.*conj(FFTX))
            MX = MX*2;                              % Multiply by 2 to take into account the fact that we threw out second half of FFTX above
            MX(1) = MX(1)/2;                        % Account for endpoint uniqueness
            MX(length(MX)) = MX(length(MX))/2;      % We know NFFT is even
            MX = MX/length(signal);                 % Scale the FFT so that it is not a function of the length of x.
            MX = MX.^2;
            MX_all(:, n) = sqrt(MX);
        end
    end
    f = f';
    mean_MX = mean(MX_all, 2);
    sem_MX = std(MX_all, 0, 2)/sqrt(size(MX_all, 2));
    figure;
    patch([f; f(end:-1:1)], [mean_MX+sem_MX; mean_MX(end:-1:1)-sem_MX(end:-1:1)], [0.8 0.8 0.8], 'EdgeColor', [0.8 0.8 0.8]);
    hold on;
    plot(f, mean_MX, '-k');
    xlim([0 FrameRate/2]);
    xlabel('Frequency (Hz)');
    ylabel('Amplitude (mm)');
    title('Z');
    box off;
catch
    errordlg('Cannot perform Fourier transform for Z positions');
end

t = (-trange:1:trange)'/FrameRate;
time_range = [t(1) t(end)];
if ~app.NoLightButton.Value
    for i = 1:ninhibition
        % trajectories aligned to bite time
        figure;
        plot_tj_multitrial(time_range, t, zbite_NI(:, 1:nbite_NI(i), i), t, zbite_I(:, 1:nbite_I(i), i), 0, NaN, 'Z (mm)', ['Inhibition #' num2str(i)]);
        figure;
        plot_tj_multitrial(time_range, t, zspeedbite_NI(:, 1:nbite_NI(i), i), t, zspeedbite_I(:, 1:nbite_I(i), i), 0, NaN, 'Speed (mm/s)', ['Z (Inhibition #' num2str(i) ')']);
        figure;
        plot_tj_multitrial(time_range, t, xyzspeedbite_NI(:, 1:nbite_NI(i), i), t, xyzspeedbite_I(:, 1:nbite_I(i), i), 0, NaN, 'Speed (mm/s)', ['XYZ (Inhibition #' num2str(i) ')']);
        figure;
        plot_tj_multitrial(time_range, t, zaccelerationbite_NI(:, 1:nbite_NI(i), i), t, zaccelerationbite_I(:, 1:nbite_I(i), i), 0, NaN, 'Absolute Acceleration (mm/s^{2})', ['Z (Inhibition #' num2str(i) ')']);
        figure;
        plot_tj_multitrial(time_range, t, xyzaccelerationbite_NI(:, 1:nbite_NI(i), i), t, xyzaccelerationbite_I(:, 1:nbite_I(i), i), 0, NaN, 'Absolute Acceleration (mm/s^{2})', ['XYZ (Inhibition #' num2str(i) ')']);
        figure;
        plot_tj_multitrial(time_range, t, abs(zphasebite_NI(:, 1:nbite_NI(i), i)), t, abs(zphasebite_I(:, 1:nbite_I(i), i)), 0, NaN, ['Absolute Phase (' char(176) ')'], ['Inhibition #' num2str(i)]);
        
        % histograms
        if nbite_I(i) ~= 0
            histogram_multitrial(zlight_all_NI{i}, zbite_NI(trange+1, 1:nbite_NI(i), i), zlight_all_I{i}, zbite_I(trange+1, 1:nbite_I(i), i), 'Z (mm)', ['Z (Inhibition #' num2str(i) ')']);
            histogram_multitrial(zspeedlight_all_NI{i}, zspeedbite_NI(trange+1, 1:nbite_NI(i), i), zspeedlight_all_I{i}, zspeedbite_I(trange+1, 1:nbite_I(i), i), 'Speed (mm/s)', ['Z (Inhibition #' num2str(i) ')']);
            histogram_multitrial(xyzspeedlight_all_NI{i}, xyzspeedbite_NI(trange+1, 1:nbite_NI(i), i), xyzspeedlight_all_I{i}, xyzspeedbite_I(trange+1, 1:nbite_I(i), i), 'Speed (mm/s)', ['XYZ (Inhibition #' num2str(i) ')']);
            histogram_multitrial(zaccelerationlight_all_NI{i}, zaccelerationbite_NI(trange+1, 1:nbite_NI(i), i), zaccelerationlight_all_I{i}, zaccelerationbite_I(trange+1, 1:nbite_I(i), i), 'Absolute Acceleration (mm/s^{2})', ['Z (Inhibition #' num2str(i) ')']);
            histogram_multitrial(xyzaccelerationlight_all_NI{i}, xyzaccelerationbite_NI(trange+1, 1:nbite_NI(i), i), xyzaccelerationlight_all_I{i}, xyzaccelerationbite_I(trange+1, 1:nbite_I(i), i), 'Absolute Acceleration (mm/s^{2})', ['XYZ (Inhibition #' num2str(i) ')']);
            
            histogram_multitrial(zphasebite_NI(trange+1, 1:nbite_NI(i), i), [], zphasebite_I(trange+1, 1:nbite_I(i), i), [], ['Phase (' char(176) ')'], ['Phase@bite (Inhibition #' num2str(i) ')']);
        else
            histogram_multitrial(zlight_all_NI{i}, zbite_NI(trange+1, 1:nbite_NI(i), i), zlight_all_I{i}, [], 'Z (mm)', ['Z (Inhibition #' num2str(i) ')']);
            histogram_multitrial(zspeedlight_all_NI{i}, zspeedbite_NI(trange+1, 1:nbite_NI(i), i), zspeedlight_all_I{i}, [], 'Speed (mm/s)', ['Z (Inhibition #' num2str(i) ')']);
            histogram_multitrial(xyzspeedlight_all_NI{i}, xyzspeedbite_NI(trange+1, 1:nbite_NI(i), i), xyzspeedlight_all_I{i}, [], 'Speed (mm/s)', ['XYZ (Inhibition #' num2str(i) ')']);
            histogram_multitrial(zaccelerationlight_all_NI{i}, zaccelerationbite_NI(trange+1, 1:nbite_NI(i), i), zaccelerationlight_all_I{i}, [], 'Absolute Acceleration (mm/s^{2})', ['Z (Inhibition #' num2str(i) ')']);
            histogram_multitrial(xyzaccelerationlight_all_NI{i}, xyzaccelerationbite_NI(trange+1, 1:nbite_NI(i), i), xyzaccelerationlight_all_I{i}, [], 'Absolute Acceleration (mm/s^{2})', ['XYZ (Inhibition #' num2str(i) ')']);
        end
        
        histogram_multitrial(log(zlight_var_NI{i}), [], log(zlight_var_I{i}), [], 'Log Z Variance (mm^{2})', ['Z (Inhibition #' num2str(i) ')']);
        histogram_multitrial(log(zspeedlight_var_NI{i}), [], log(zspeedlight_var_I{i}), [], 'Log Speed Variance (mm^{2}/s^{2})', ['Z (Inhibition #' num2str(i) ')']);
        histogram_multitrial(log(xyzspeedlight_var_NI{i}), [], log(xyzspeedlight_var_I{i}), [], 'Log Speed Variance (mm^{2}/s^{2})', ['XYZ (Inhibition #' num2str(i) ')']);
        histogram_multitrial(log(zaccelerationlight_var_NI{i}), [], log(zaccelerationlight_var_I{i}), [], 'Log Absolute Acceleration Variance (mm^{2}/s^{4})', ['Z (Inhibition #' num2str(i) ')']);
        histogram_multitrial(log(xyzaccelerationlight_var_NI{i}), [], log(xyzaccelerationlight_var_I{i}), [], 'Log Absolute Acceleration Variance (mm^{2}/s^{4})', ['XYZ (Inhibition #' num2str(i) ')']);
        
        result.zlight_mean_NI(i) = mean(zlight_mean_NI{i}, 'omitnan');
        result.zspeedlight_mean_NI(i) = mean(zspeedlight_mean_NI{i}, 'omitnan');
        result.xyzspeedlight_mean_NI(i) = mean(xyzspeedlight_mean_NI{i}, 'omitnan');
        result.zaccelerationlight_mean_NI(i) = mean(zaccelerationlight_mean_NI{i}, 'omitnan');
        result.xyzaccelerationlight_mean_NI(i) = mean(xyzaccelerationlight_mean_NI{i}, 'omitnan');
        
        result.zlight_var_NI(i) = mean(zlight_var_NI{i}, 'omitnan');
        result.zspeedlight_var_NI(i) = mean(zspeedlight_var_NI{i}, 'omitnan');
        result.xyzspeedlight_var_NI(i) = mean(xyzspeedlight_var_NI{i}, 'omitnan');
        result.zaccelerationlight_var_NI(i) = mean(zaccelerationlight_var_NI{i}, 'omitnan');
        result.xyzaccelerationlight_var_NI(i) = mean(xyzaccelerationlight_var_NI{i}, 'omitnan');
        
        result.zlight_mean_I(i) = mean(zlight_mean_I{i}, 'omitnan');
        result.zspeedlight_mean_I(i) = mean(zspeedlight_mean_I{i}, 'omitnan');
        result.xyzspeedlight_mean_I(i) = mean(xyzspeedlight_mean_I{i}, 'omitnan');
        result.zaccelerationlight_mean_I(i) = mean(zaccelerationlight_mean_I{i}, 'omitnan');
        result.xyzaccelerationlight_mean_I(i) = mean(xyzaccelerationlight_mean_I{i}, 'omitnan');
        
        result.zlight_var_I(i) = mean(zlight_var_I{i}, 'omitnan');
        result.zspeedlight_var_I(i) = mean(zspeedlight_var_I{i}, 'omitnan');
        result.xyzspeedlight_var_I(i) = mean(xyzspeedlight_var_I{i}, 'omitnan');
        result.zaccelerationlight_var_I(i) = mean(zaccelerationlight_var_I{i}, 'omitnan');
        result.xyzaccelerationlight_var_I(i) = mean(xyzaccelerationlight_var_I{i}, 'omitnan');
        
        result.zbite_NI{i} = zbite_NI(:, 1:nbite_NI(i), i);
        result.zspeedbite_NI{i} = zspeedbite_NI(:, 1:nbite_NI(i), i);
        result.xyzspeedbite_NI{i} = xyzspeedbite_NI(:, 1:nbite_NI(i), i);
        result.zaccelerationbite_NI{i} = zaccelerationbite_NI(:, 1:nbite_NI(i), i);
        result.xyzaccelerationbite_NI{i} = xyzaccelerationbite_NI(:, 1:nbite_NI(i), i);
        
        result.zbite_I{i} = zbite_I(:, 1:nbite_I(i), i);
        result.zspeedbite_I{i} = zspeedbite_I(:, 1:nbite_I(i), i);
        result.xyzspeedbite_I{i} = xyzspeedbite_I(:, 1:nbite_I(i), i);
        result.zaccelerationbite_I{i} = zaccelerationbite_I(:, 1:nbite_I(i), i);
        result.xyzaccelerationbite_I{i} = xyzaccelerationbite_I(:, 1:nbite_I(i), i);
        
        result.zphasebite_NI{i} = zphasebite_NI(:, 1:nbite_NI(i), i);
        result.zphasebite_I{i} = zphasebite_I(:, 1:nbite_I(i), i);
    end
else
    for i = 1:ninhibition
        % trajectories aligned to bite time
        figure;
        plot_tj_multitrial(time_range, t, zbite_NI(:, 1:nbite_NI(i), i), [], [], 0, NaN, 'Z (mm)', ['Inhibition #' num2str(i)]);
        figure;
        plot_tj_multitrial(time_range, t, zspeedbite_NI(:, 1:nbite_NI(i), i), [], [], 0, NaN, 'Speed (mm/s)', ['Z (Inhibition #' num2str(i) ')']);
        figure;
        plot_tj_multitrial(time_range, t, xyzspeedbite_NI(:, 1:nbite_NI(i), i), [], [], 0, NaN, 'Speed (mm/s)', ['XYZ (Inhibition #' num2str(i) ')']);
        figure;
        plot_tj_multitrial(time_range, t, zaccelerationbite_NI(:, 1:nbite_NI(i), i), [], [], 0, NaN, 'Absolute Acceleration (mm/s^{2})', ['Z (Inhibition #' num2str(i) ')']);
        figure;
        plot_tj_multitrial(time_range, t, xyzaccelerationbite_NI(:, 1:nbite_NI(i), i), [], [], 0, NaN, 'Absolute Acceleration (mm/s^{2})', ['XYZ (Inhibition #' num2str(i) ')']);
        if app.sDelay4sOnButton.Value
            figure;
            plot_tj_multitrial(time_range, t, abs(zphasebite_NI(:, 1:nbite_NI(i), i)), [], [], 0, NaN, ['Absolute Phase (' char(176) ')'], ['Inhibition #' num2str(i)]);
        end
        
        % histograms
        histogram_multitrial(zlight_all_NI{i}, zbite_NI(trange+1, 1:nbite_NI(i), i), [], [], 'Z (mm)', ['Z (Inhibition #' num2str(i) ')']);
        histogram_multitrial(zspeedlight_all_NI{i}, zspeedbite_NI(trange+1, 1:nbite_NI(i), i), [], [], 'Speed (mm/s)', ['Z (Inhibition #' num2str(i) ')']);
        histogram_multitrial(xyzspeedlight_all_NI{i}, xyzspeedbite_NI(trange+1, 1:nbite_NI(i), i), [], [], 'Speed (mm/s)', ['XYZ (Inhibition #' num2str(i) ')']);
        histogram_multitrial(zaccelerationlight_all_NI{i}, zaccelerationbite_NI(trange+1, 1:nbite_NI(i), i), [], [], 'Absolute Acceleration (mm/s^{2})', ['Z (Inhibition #' num2str(i) ')']);
        histogram_multitrial(xyzaccelerationlight_all_NI{i}, xyzaccelerationbite_NI(trange+1, 1:nbite_NI(i), i), [], [], 'Absolute Acceleration (mm/s^{2})', ['XYZ (Inhibition #' num2str(i) ')']);
        
        histogram_multitrial(log(zlight_var_NI{i}), [], [], [], 'Log Z Variance (mm^{2})', ['Z (Inhibition #' num2str(i) ')']);
        histogram_multitrial(log(zspeedlight_var_NI{i}), [], [], [], 'Log Speed Variance (mm^{2}/s^{2})', ['Z (Inhibition #' num2str(i) ')']);
        histogram_multitrial(log(xyzspeedlight_var_NI{i}), [], [], [], 'Log Speed Variance (mm^{2}/s^{2})', ['XYZ (Inhibition #' num2str(i) ')']);
        histogram_multitrial(log(zaccelerationlight_var_NI{i}), [], [], [], 'Log Absolute Acceleration Variance (mm^{2}/s^{4})', ['Z (Inhibition #' num2str(i) ')']);
        histogram_multitrial(log(xyzaccelerationlight_var_NI{i}), [], [], [], 'Log Absolute Acceleration Variance (mm^{2}/s^{4})', ['XYZ (Inhibition #' num2str(i) ')']);
        
        result.zlight_mean_NI(i) = mean(zlight_mean_NI{i}, 'omitnan');
        result.zspeedlight_mean_NI(i) = mean(zspeedlight_mean_NI{i}, 'omitnan');
        result.xyzspeedlight_mean_NI(i) = mean(xyzspeedlight_mean_NI{i}, 'omitnan');
        result.zaccelerationlight_mean_NI(i) = mean(zaccelerationlight_mean_NI{i}, 'omitnan');
        result.xyzaccelerationlight_mean_NI(i) = mean(xyzaccelerationlight_mean_NI{i}, 'omitnan');
        
        result.zlight_var_NI(i) = mean(zlight_var_NI{i}, 'omitnan');
        result.zspeedlight_var_NI(i) = mean(zspeedlight_var_NI{i}, 'omitnan');
        result.xyzspeedlight_var_NI(i) = mean(xyzspeedlight_var_NI{i}, 'omitnan');
        result.zaccelerationlight_var_NI(i) = mean(zaccelerationlight_var_NI{i}, 'omitnan');
        result.xyzaccelerationlight_var_NI(i) = mean(xyzaccelerationlight_var_NI{i}, 'omitnan');
        
        result.zbite_NI{i} = zbite_NI(:, 1:nbite_NI(i), i);
        result.zspeedbite_NI{i} = zspeedbite_NI(:, 1:nbite_NI(i), i);
        result.xyzspeedbite_NI{i} = xyzspeedbite_NI(:, 1:nbite_NI(i), i);
        result.zaccelerationbite_NI{i} = zaccelerationbite_NI(:, 1:nbite_NI(i), i);
        result.xyzaccelerationbite_NI{i} = xyzaccelerationbite_NI(:, 1:nbite_NI(i), i);
        
        result.zbite_NI_1p{i} = zbite_NI(trange+1, 1:nbite_NI(i), i);
        
        if app.sDelay4sOnButton.Value
            result.zphasebite_NI{i} = zphasebite_NI(:, 1:nbite_NI(i), i);
        end
    end
end

curpwd = pwd;
try
    cd(path);
end
[file, path] = uiputfile('AverageTrajectory.mat', 'Save result');
if file ~= 0
    save([path file], 'result');
end
cd(curpwd);