%*---------------------------------------------------------------------*
% Huang Lab
% Cold Spring Harbor Laboratory
% Author : Xu An, June 2019
% xan@cshl.edu
% Version: 1.0
%*---------------------------------------------------------------------*
%%
%load raw .bin file into workspace
%% batch processing
optsession = 1;
savemat = 1;

pathname = {'I:\Free-moving Feeding Data\20221010_AngelHair_Fezf2CreER;PlexinD1flp_Gcamp-R1',...
    'I:\Free-moving Feeding Data\20221010_AngelHair_Fezf2CreER;PlexinD1flp_Gcamp-R2',...
    'I:\Free-moving Feeding Data\20221013_AngelHair_Fezf2CreER;PlexinD1flp_Gcamp-R1',...
    'I:\Free-moving Feeding Data\20221013_AngelHair_Fezf2CreER;PlexinD1flp_Gcamp-R2',...
    'I:\Free-moving Feeding Data\20221014_AngelHair_Fezf2CreER;PlexinD1flp_Gcamp-R1',...
    'I:\Free-moving Feeding Data\20221014_AngelHair_Fezf2CreER;PlexinD1flp_Gcamp-R2'};

% pathname = {'I:\Free-moving Feeding Data\20210211_AngelHair_Fezf2CreER_Gcamp8'};
    
filename = '\Session001.bin';

hf = figure('Color', [1 1 1]);
workbar(0, 'Computing Ongoing...', 'Progress');
for ii = 1:numel(pathname)
    clf(hf);
    set(hf, 'Name', pathname{ii});
    clear A;
    fileID = fopen([pathname{ii} filename]);
    if optsession
        A = fread(fileID, [4 inf], 'double');
    else
        A = fread(fileID, [3 inf], 'double');
    end
    fclose(fileID);
    trial_mark = diff(A(2, :));
    trial_starts = find(trial_mark == 1)+1;
    trial_ends = find(trial_mark == -1);
    if optsession
        laser_mark = diff(A(3, :));
        laser_starts = find(laser_mark == 1)+1;
        laser_ends = find(laser_mark == -1);
    end
    
    totaltrial = numel(trial_starts);
    for i = 1:totaltrial
        audiodata = [];
        timestamps = A(1, trial_starts(i):trial_ends(i));
        if optsession
            laser_starttime = [];
            laser_endtime = [];
            laser_index = find(laser_starts >= trial_starts(i) & laser_starts <= trial_ends(i));
            if ~isempty(laser_index)
                laser_starttime = A(1, laser_starts(laser_index));
                laser_endtime = A(1, laser_ends(laser_index));
            end
            laser_starttime = laser_starttime-timestamps(1);
            laser_endtime = laser_endtime-timestamps(1);
            timestamps = timestamps-timestamps(1);
            audiosignal = A(4, trial_starts(i):trial_ends(i));
        else
            timestamps = timestamps-timestamps(1);
            audiosignal = A(3, trial_starts(i):trial_ends(i));
        end
        
        if i <= 48
            haxis = subplot(8, 6, i, 'Parent', hf);
        else
            figure;
            haxis = axes;
        end
        if optsession
            if ~isempty(laser_index)
                for j = 1:numel(laser_index)
                    patch(haxis, [laser_starttime(j) laser_starttime(j) laser_endtime(j) laser_endtime(j)], [min(audiosignal) max(audiosignal) max(audiosignal) min(audiosignal)], [0 1 0], 'FaceAlpha', 1, 'EdgeColor', 'none');
                end
                hold(haxis, 'on');
            end
        end
        plot(haxis, timestamps, audiosignal, '-k');
        set(haxis, 'xLim', [0 timestamps(end)]);
        set(haxis, 'yLim', [min(audiosignal) max(audiosignal)]);
        xlabel(haxis, 'Time (s)');
        ylabel(haxis, 'Voltage (V)');
        box(haxis, 'off');
        title(haxis, ['Trial ' num2str(i)]);
        drawnow;
        
        if savemat
            audiodata.timestamps = timestamps;
            if optsession
                audiodata.laser_starttime = laser_starttime;
                audiodata.laser_endtime = laser_endtime;
            end
            audiodata.signal = audiosignal;
            save([pathname{ii} '\Trial' num2str(i, '%03d') '.mat'], 'audiodata');
        end
    end
    workbar(ii/numel(pathname), [num2str(ii) '/' num2str(numel(pathname))], 'Progress');
end
msgbox('Done !');

%% interactively
optsession = 1;
clear A;
curpwd = pwd;
try
    cd(pathname);
end
[filename, pathname] = uigetfile({'*.bin'}, 'Select audio data');
if filename ~= 0
    fileID = fopen([pathname filename]);
    if optsession
        A = fread(fileID, [4 inf], 'double');
    else
        A = fread(fileID, [3 inf], 'double');
    end
    fclose(fileID);
    trial_mark = diff(A(2, :));
    trial_starts = find(trial_mark == 1)+1;
    trial_ends = find(trial_mark == -1);
    if optsession
        laser_mark = diff(A(3, :));
        laser_starts = find(laser_mark == 1)+1;
        laser_ends = find(laser_mark == -1);
    end
end
cd(curpwd);
msgbox('Done !');

%%
%plot out all trials and save them into .mat files
savemat = 1;
totaltrial = numel(trial_starts);
workbar(0, 'Computing Ongoing...', 'Progress');
hf = figure('Color', [1 1 1]);
for i = 1:totaltrial
    audiodata = [];
    timestamps = A(1, trial_starts(i):trial_ends(i));
    if optsession
        laser_starttime = [];
        laser_endtime = [];
        laser_index = find(laser_starts >= trial_starts(i) & laser_starts <= trial_ends(i));
        if ~isempty(laser_index)
            laser_starttime = A(1, laser_starts(laser_index));
            laser_endtime = A(1, laser_ends(laser_index));
        end
        laser_starttime = laser_starttime-timestamps(1);
        laser_endtime = laser_endtime-timestamps(1);
        timestamps = timestamps-timestamps(1);
        audiosignal = A(4, trial_starts(i):trial_ends(i));
    else
        timestamps = timestamps-timestamps(1);
        audiosignal = A(3, trial_starts(i):trial_ends(i));
    end
    
    haxis = subplot(8, 6, i, 'Parent', hf);
    if optsession
        if ~isempty(laser_index)
            for j = 1:numel(laser_index)
                patch(haxis, [laser_starttime(j) laser_starttime(j) laser_endtime(j) laser_endtime(j)], [min(audiosignal) max(audiosignal) max(audiosignal) min(audiosignal)], [0 1 0], 'FaceAlpha', 1, 'EdgeColor', 'none');
            end
            hold(haxis, 'on');
        end
    end
    plot(haxis, timestamps, audiosignal, '-k');
    set(haxis, 'xLim', [0 timestamps(end)], 'ButtonDownFcn', @extract_figure);
    set(haxis, 'yLim', [min(audiosignal) max(audiosignal)]);
    xlabel(haxis, 'Time (s)');
    ylabel(haxis, 'Voltage (V)');
    box(haxis, 'off');
    title(haxis, ['Trial ' num2str(i)]);
    
    if savemat
        audiodata.timestamps = timestamps;
        if optsession
            audiodata.laser_starttime = laser_starttime;
            audiodata.laser_endtime = laser_endtime;
        end
        audiodata.signal = audiosignal;
        save([pathname 'Trial' num2str(i, '%03d') '.mat'], 'audiodata');
    end
    workbar(i/totaltrial, ['Trial' num2str(i) '/' num2str(totaltrial)], 'Progress');
end

%%
%play the sound from the last trial
SampleRate = 96000;
soundsc(audiosignal, SampleRate);

%%
%low pass filter audio data of single trial and save the result as a .wav file to play
%outside matlab
SampleRate = 96000;
curpwd = pwd;
try
    cd(pathname);
end
[filename, pathname] = uigetfile({'*.mat'}, 'Select audio data');
if filename ~= 0
    load([pathname filename]);
else
    cd(curpwd);
    return;
end

[filename, pathname] = uiputfile({[filename(1:end-4) '.wav']}, 'Save audio data as');
if filename ~= 0
    freq = 8000;
    [b, a] = butter(6, freq/(SampleRate/2), 'low');
    signal_filtered = filtfilt(b, a, audiodata.signal);
    audiowrite([pathname filename], rescale(signal_filtered, -1, 1), SampleRate, 'BitsPerSample', 16)
end
cd(curpwd);

%%
%write audio data into the video file
curpwd = pwd;
try
    cd(pathname);
end
[filename, pathname] = uigetfile({'*.avi'}, 'Select video data');
if filename ~= 0
    cd(pathname);
    SampleRate = 96000;
    freq = 8000;
    FrameRate_scale = 0.5;
    sqsize = 90;
    
    videoCSVname = [filename(1:3) 'gpio' filename(10:end-4) '.csv'];
    videoCSVdata = load(videoCSVname);
    framecount = videoCSVdata(:, 2);
    framecount = framecount-framecount(1)+1;
    
    cd ..
    load(['Trial' num2str(str2double(filename(10:end-4)), '%03d') '.mat']);
    [b, a] = butter(6, freq/(SampleRate/2), 'low');
    signal_filtered = filtfilt(b, a, audiodata.signal);
    signal_filtered = rescale(signal_filtered);
    
    vidObj = VideoReader([pathname filename]);
    FrameRate = vidObj.FrameRate;
    time_interval = 1/FrameRate;
    FrameRate_save = FrameRate*FrameRate_scale;
    videoFReader = vision.VideoFileReader([pathname filename]);
    
    videoFWriter = vision.VideoFileWriter(filename, 'FrameRate', FrameRate_save, 'AudioInputPort', 1, 'VideoCompressor', 'MJPEG Compressor');
    n = 1;
    while ~isDone(videoFReader)
        videoFrame = step(videoFReader);
        if videoCSVdata(n, 1)
            videoFrame = insertShape(videoFrame, 'FilledRectangle', [0 0 sqsize sqsize], 'Color', 'w', 'Opacity', 1);
        end
        start_index = find(audiodata.timestamps >= (framecount(n)-1)*time_interval, 1);
        end_index = start_index+round(SampleRate*time_interval)-1;
        if end_index > numel(signal_filtered)
            break;
        end
        audioFrame = signal_filtered(start_index:end_index);
        step(videoFWriter, videoFrame, audioFrame');
        n = n+1;
    end
    release(videoFReader);
    release(videoFWriter);
end
cd(curpwd);
msgbox('Done !');