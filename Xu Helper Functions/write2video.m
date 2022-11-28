function write2video(app, FrameRate, csvdata, readerobj, num_track, duration, label_table, label_trace, sqsize, dotsize, labelcolor, alphavalue, fpdata)
brightness = app.BrightnessEditField.Value;
videoFReader = vision.VideoFileReader([readerobj.Path '\' readerobj.Name]);
videoFWriter = vision.VideoFileWriter([readerobj.Path '\Short_' readerobj.Name(1:end-4) '.mp4'],...
    'FrameRate', FrameRate, 'FileFormat', 'MPEG4', 'Quality', 100, 'VideoCompressor', 'MJPEG Compressor');
n = 1;
frames2stop = duration(2)*readerobj.FrameRate;
nlabels = size(label_table, 1);
ltrace = NaN(nlabels, 2, label_trace);

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
while ~isDone(videoFReader) && n <= frames2stop
    videoFrame = step(videoFReader);
    if ~isempty(num_track)
        for i = 1:nlabels
            if label_table(i, 2) && num_track(n, i*3+1) > label_table(i, 1)
                ltrace(i, 1, end) = num_track(n, i*3-1);
                ltrace(i, 2, end) = num_track(n, i*3);
            else
                ltrace(i, 1, end) = NaN;
                ltrace(i, 2, end) = NaN;
            end
        end
    end
    if n >= duration(1)*readerobj.FrameRate
        if brightness ~= 1
            videoFrame = rgb2hsv(videoFrame);
            intensity = videoFrame(:, :, 3);
            intensity = mat2gray(intensity, [0 brightness]);
            videoFrame(:, :, 3) = intensity;
            videoFrame = hsv2rgb(videoFrame);
        end
        if csvdata(n, 1) && app.OptSessionCheckBox.Value
            videoFrame = insertShape(videoFrame, 'FilledRectangle', [0 0 sqsize sqsize], 'Color', 'w', 'Opacity', 1);
        end
        if ~isempty(num_track)
            for i = 1:nlabels
                for j = 1:label_trace
                    if ~isnan(ltrace(i, 1, j))
                        videoFrame = insertShape(videoFrame, 'FilledCircle', [ltrace(i, 1, j) ltrace(i, 2, j) dotsize], 'Color',...
                            labelcolor(round((i-1)*255/(nlabels-1))+1, :)/(2^8), 'Opacity', j/label_trace*alphavalue);
                    end
                end
            end
        end
        
        if ~isempty(fpdata) && app.FPCheckBox.Value
            if ~ha_im.UserData
                him = imshow(videoFrame, 'Parent', ha_im);
                ha_im.UserData = 1;
            else
                him.CData = videoFrame;
            end
            hpatch.XData = [(n+1)/readerobj.FrameRate (n+1)/readerobj.FrameRate duration(2) duration(2)];
            videoFrame = getframe(hf);
            videoFrame = videoFrame.cdata;
        end
        step(videoFWriter, videoFrame);
    end
    ltrace = circshift(ltrace, -1, 3);
    workbar(n/floor(frames2stop), [num2str(n) '/' num2str(floor(frames2stop))], 'Progress');
    n = n+1;
end
release(videoFReader);
release(videoFWriter);
if ~isempty(fpdata) && app.FPCheckBox.Value
    close(hf);
end
end