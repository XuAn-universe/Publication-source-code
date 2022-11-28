function writesound2video2(app, FrameRate, SampleRate, audiodata, csvdata, readerobj, num_track, signal_filtered, duration, label_table, label_trace, sqsize, dotsize, labelcolor, alphavalue, fpdata)
brightness = app.BrightnessEditField.Value;
framecount = csvdata{2}(:, 2);
framecount = framecount-framecount(1)+1;
time_interval = 1/readerobj{2}.FrameRate;
sideID = 1;
videoFReader1 = vision.VideoFileReader([readerobj{sideID}.Path '\' readerobj{sideID}.Name]);
videoFReader2 = vision.VideoFileReader([readerobj{2}.Path '\' readerobj{2}.Name]);
videoFWriter = vision.VideoFileWriter([readerobj{2}.Path '\Sound_' readerobj{sideID}.Name(1:end-4) readerobj{2}.Name],...
    'FrameRate', FrameRate, 'AudioInputPort', 1, 'VideoCompressor', 'MJPEG Compressor'); % None (uncompressed)
n = 1;
frames2stop = duration(2)*readerobj{2}.FrameRate;
nlabels = size(label_table, 1);
ltrace1 = NaN(nlabels, 2, label_trace);
ltrace2 = NaN(nlabels, 2, label_trace);

% get photometry data
if ~isempty(fpdata) && app.FPCheckBox.Value
    value = app.TrialsListBox.Value;
    % get bite events
    Bite_events = [];
    try
        audiolocation = app.EditField.Value(1:end-7);
        temp = load([audiolocation '\Detected_Bite_Events.mat']);
        Bite_events = temp.Audio_analysis;
    end
    bite_timestamps = [];
    if ~isempty(Bite_events)
        bite_timestamps = Bite_events(value).time_bites;
        if ~isempty(bite_timestamps)
%             bite_amplitudes = Bite_events(value).amplitude_bites/(max(Bite_events(value).amplitude_bites));
            bite_amplitudes = Bite_events(value).amplitude_bites./Bite_events(value).amplitude_bites;
        end
    end
    
    LabelledEvents = [];
    try
        temp = load([app.EditField.Value '\LabelledEvents' num2str(value) '.mat']);
        LabelledEvents = temp.LabelledEvents;
    end

    nchannel = numel(fpdata);
    fpdata_t = fpdata{1}(:, 1);
    for i = 1:nchannel
        fpdata_zsignal(:, i) = fpdata{i}(:, 2);
    end
    [yl(1), yl(2)] = bounds(fpdata_zsignal(fpdata_t >= duration(1) & fpdata_t <= duration(2), :), 'all');
    zscore2y = 6/(yl(2)-yl(1));
    zero_position = zscore2y*(-yl(1));
    
    vWidth = round(readerobj.Width/2);
    vHeight = round(readerobj.Height/2);
    width_scale = 5; % must be >= 1
    height_scale = 1;
    plot_width = round(vWidth*width_scale);
    plot_height = round(vHeight*height_scale);
    hf = figure('Units', 'pixels', 'Position', [0, 0, plot_width, vHeight+plot_height], 'Color', [1 1 1], 'Visible', 'off');
    ha_im = axes('Units', 'pixels', 'Position', [round((plot_width-vWidth)/2) plot_height vWidth, vHeight], 'UserData', 0);
    ha_plot = axes('Units', 'pixels', 'OuterPosition', [0 0 plot_width, plot_height]);
    hold(ha_plot, 'on');
    
    if ~isempty(LabelledEvents)
        plot_events(ha_plot, 1, 1, LabelledEvents);
    end
    
    if ~isempty(bite_timestamps)
        for i = 1:numel(bite_timestamps)
            line(ha_plot, [bite_timestamps(i) bite_timestamps(i)], [0.5-bite_amplitudes(i)/2 0.5+bite_amplitudes(i)/2],...
                'Color', [1 0 1], 'LineStyle', '-', 'LineWidth', 1);
        end
    end
    
    for i = 1:nchannel
        hp(i) = plot(ha_plot, fpdata_t, zero_position+fpdata_zsignal(:, i)*zscore2y, 'Color', [1-1/nchannel*i 1-1/nchannel*i 1-1/nchannel*i], 'LineWidth', 1);
    end
    line(ha_plot, [fpdata_t(1), fpdata_t(end)], [zero_position zero_position], 'Color', [0 0 0], 'LineStyle', '--', 'LineWidth', 0.5);
    
    hpatch = patch(ha_plot, 'XData', [duration(1) duration(1) duration(2) duration(2)], 'YData', [0+0.03 6 6 0+0.03], 'FaceColor', [1 1 1], 'EdgeColor', 'none');

    lgdtext{1} = 'aCFA'; %
    lgdtext{2} = 'RFO'; %
    xlim(ha_plot, duration);
    ylim(ha_plot, [0 6]);
    set(ha_plot, 'YTick', [0 zero_position 6], 'YTickLabel', {num2str(yl(1), '%.1f'), '0', num2str(yl(2), '%.1f')}, 'TickLength', [0 0], 'FontSize', 12);
    xlabel(ha_plot, 'Time (s)');
    ylabel(ha_plot, 'Z-score');
    legend(hp, lgdtext, 'Location', 'northeastoutside', 'FontSize', 12);
    legend('boxoff');
    box(ha_plot, 'off');
end

workbar(0, 'Computing Ongoing...', 'Progress');
while ~isDone(videoFReader2) && n <= frames2stop
    videoFrame1 = step(videoFReader1);
    videoFrame2 = step(videoFReader2);
    if ~isempty(num_track{sideID})
        for i = 1:nlabels
            if label_table(i, 2) && num_track{sideID}(n, i*3+1) > label_table(i, 1)
                ltrace1(i, 1, end) = num_track{sideID}(n, i*3-1);
                ltrace1(i, 2, end) = num_track{sideID}(n, i*3);
            else
                ltrace1(i, 1, end) = NaN;
                ltrace1(i, 2, end) = NaN;
            end
        end
    end
    if ~isempty(num_track{2})
        for i = 1:nlabels
            if label_table(i, 2) && num_track{2}(n, i*3+1) > label_table(i, 1)
                ltrace2(i, 1, end) = num_track{2}(n, i*3-1);
                ltrace2(i, 2, end) = num_track{2}(n, i*3);
            else
                ltrace2(i, 1, end) = NaN;
                ltrace2(i, 2, end) = NaN;
            end
        end
    end
    if n >= duration(1)*readerobj{2}.FrameRate
        if brightness ~= 1
            videoFrame1 = rgb2hsv(videoFrame1);
            intensity1 = videoFrame1(:, :, 3);
            intensity1 = mat2gray(intensity1, [0 brightness]);
            videoFrame1(:, :, 3) = intensity1;
            videoFrame1 = hsv2rgb(videoFrame1);
            
            videoFrame2 = rgb2hsv(videoFrame2);
            intensity2 = videoFrame2(:, :, 3);
            intensity2 = mat2gray(intensity2, [0 brightness]);
            videoFrame2(:, :, 3) = intensity2;
            videoFrame2 = hsv2rgb(videoFrame2);
        end
        if ~app.sOnTwiceButton.Value
            if csvdata{2}(n, 1) && app.OptSessionCheckBox.Value
                videoFrame1 = insertShape(videoFrame1, 'FilledRectangle', [0 0 sqsize sqsize], 'Color', 'w', 'Opacity', 1);
                
                videoFrame2 = insertShape(videoFrame2, 'FilledRectangle', [0 0 sqsize sqsize], 'Color', 'w', 'Opacity', 1);
            end
        else
            if ~isempty(audiodata.laser_starttime) && app.OptSessionCheckBox.Value
                if ((n/readerobj{2}.FrameRate >= audiodata.laser_starttime(1) && n/readerobj{2}.FrameRate <= audiodata.laser_endtime(1))...
                        || (n/readerobj{2}.FrameRate >= audiodata.laser_starttime(2) && n/readerobj{2}.FrameRate <= audiodata.laser_endtime(2)))
                    videoFrame1 = insertShape(videoFrame1, 'FilledRectangle', [0 0 sqsize sqsize], 'Color', 'w', 'Opacity', 1);
                    
                    videoFrame2 = insertShape(videoFrame2, 'FilledRectangle', [0 0 sqsize sqsize], 'Color', 'w', 'Opacity', 1);
                end
            end
        end
        if ~isempty(num_track{sideID})
            for i = 1:nlabels
                for j = 1:label_trace
                    if ~isnan(ltrace1(i, 1, j))
                        videoFrame1 = insertShape(videoFrame1, 'FilledCircle', [ltrace1(i, 1, j) ltrace1(i, 2, j) dotsize], 'Color',...
                            labelcolor(round((i-1)*255/(nlabels-1))+1, :)/(2^8), 'Opacity', j/label_trace*alphavalue);
                    end
                end
            end
        end
        if ~isempty(num_track{2})
            for i = 1:nlabels
                for j = 1:label_trace
                    if ~isnan(ltrace2(i, 1, j))
                        videoFrame2 = insertShape(videoFrame2, 'FilledCircle', [ltrace2(i, 1, j) ltrace2(i, 2, j) dotsize], 'Color',...
                            labelcolor(round((i-1)*255/(nlabels-1))+1, :)/(2^8), 'Opacity', j/label_trace*alphavalue);
                    end
                end
            end
        end
        
        videoFrame = [videoFrame1 videoFrame2];
        if ~isempty(fpdata) && app.FPCheckBox.Value
            if ~ha_im.UserData
                him = imshow([videoFrame1 videoFrame2], 'Parent', ha_im);
                ha_im.UserData = 1;
            else
                him.CData = [videoFrame1 videoFrame2];
            end
            hpatch.XData = [(n+1)/readerobj.FrameRate (n+1)/readerobj.FrameRate duration(2) duration(2)];
            videoFrame = getframe(hf);
            videoFrame = videoFrame.cdata;
            
%             videoFrame = print(hf, '-RGBImage', '-r0', '-painters');
        end
        
        start_index = find(audiodata.timestamps >= (framecount(n)-1)*time_interval, 1);
        end_index = start_index+round(SampleRate*time_interval)-1;
        if end_index > numel(signal_filtered)
            break;
        end
        audioFrame = signal_filtered(start_index:end_index);
        step(videoFWriter, videoFrame, audioFrame');
    end
    ltrace1 = circshift(ltrace1, -1, 3);
    ltrace2 = circshift(ltrace2, -1, 3);
    workbar(n/floor(frames2stop), [num2str(n) '/' num2str(floor(frames2stop))], 'Progress');
    n = n+1;
end
release(videoFReader1);
release(videoFReader2);
release(videoFWriter);
if ~isempty(fpdata) && app.FPCheckBox.Value
    close(hf);
end
end