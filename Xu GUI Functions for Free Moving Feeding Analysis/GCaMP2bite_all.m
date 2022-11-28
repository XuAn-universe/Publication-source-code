function GCaMP2bite_all(app, Exp_Path)
value = app.TrialsListBox.Value;
value = sort(value);
trials = numel(value);

load([Exp_Path '\Analysis_Session.mat'], 'Video_annotation');

% get bite events
try
    audiolocation = Exp_Path(1:end-7);
    temp = load([audiolocation '\Detected_Bite_Events.mat']);
    Bite_events = temp.Audio_analysis;
catch
    errordlg('Please detect the bites first.', 'Error');
    return;
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

pre = app.presEditField.Value;
post = app.postsEditField.Value;

zscore_bite_all = [];
zscore_bite_random_all = [];
zscore_firstbite_all = [];
zscore_firstbite_random_all = [];

for i = 1:trials
    if ~Video_annotation(value(i)).Disgard
        if ~app.BiteDeviceButton.Value
            try
                temp = load([Exp_Path '\LabelledEvents' num2str(value(i)) '.mat']);
                events = temp.LabelledEvents;
            catch
                continue;
            end
        else
            events.BiteBoutStart = [];
        end
        
        bite_timestamps = Bite_events(value(i)).time_bites;
        fpdata = fpdata_all.zsignal_all(:, value(i));
        fpdata_zsignal = [];
        fpdata_t = [];
        for j = 1:nchannel
            fpdata_zsignal(:, j) = fpdata{j}(:, 2);
        end
        fpdata_t = fpdata{1}(:, 1);
        
        if ~isempty(bite_timestamps)
            if ~isempty(events.BiteBoutStart)
                bite_timestampsRnd = get_ramdomized_time_in_handling_stage4bite(events, bite_timestamps);
                [~, zscore, ~] = AlignSignal2Event(fpdata_t, fpdata_zsignal, bite_timestamps, pre, post);
                [t_aligned, zscore_random, ~] = AlignSignal2Event(fpdata_t, fpdata_zsignal, bite_timestampsRnd, pre, post);
            else
                [t_aligned, zscore, zscore_random] = AlignSignal2Event(fpdata_t, fpdata_zsignal, bite_timestamps, pre, post);
            end
            zscore_bite_all = [zscore_bite_all zscore];
            zscore_bite_random_all = [zscore_bite_random_all zscore_random];
            
            [t_aligned, zscore, zscore_random] = AlignSignal2Event(fpdata_t, fpdata_zsignal, bite_timestamps(1), pre, post);
            zscore_firstbite_all = [zscore_firstbite_all zscore];
            zscore_firstbite_random_all = [zscore_firstbite_random_all zscore_random];
        end
    end
end
figure;
for i = 1:nchannel
    subplot(1, nchannel, i);
    plot_gcamp_align2event(t_aligned, zscore_bite_random_all(:, :, i), zscore_bite_all(:, :, i), ['aligned to all bite: ' lgdtext{i}]);
end

figure;
for i = 1:nchannel
    subplot(1, nchannel, i);
    plot_gcamp_align2event(t_aligned, zscore_firstbite_random_all(:, :, i), zscore_firstbite_all(:, :, i), ['aligned to first bite: ' lgdtext{i}]);
end

for i = 1:nchannel
    figure;
    heatmap4gcamp(t_aligned, zscore_bite_all(:, :, i), ['aligned to all bite: ' lgdtext{i}]);
    
    figure;
    heatmap4gcamp(t_aligned, zscore_firstbite_all(:, :, i), ['aligned to first bite: ' lgdtext{i}]);
end

colors = [0.3010 0.7450 0.9330; 0.8500 0.3250 0.0980];
hp = zeros(1, nchannel);
figure;
subplot(1, 2, 1);
for i = 1:nchannel
    hp(i) = plot_tj_MeanSEM(t_aligned, zscore_bite_all(:, :, i), colors(i, :), colors(i, :), 'Time (s)', 'Z-score', 'aligned to all bite');
end
yl = ylim;
plot([0 0], yl, '--k', 'LineWidth', 1);
legend(hp, lgdtext, 'Location', 'NorthEast');
legend('boxoff');

subplot(1, 2, 2);
for i = 1:nchannel
    hp(i) = plot_tj_MeanSEM(t_aligned, zscore_firstbite_all(:, :, i), colors(i, :), colors(i, :), 'Time (s)', 'Z-score', 'aligned to first bite');
end
yl = ylim;
plot([0 0], yl, '--k', 'LineWidth', 1);
legend(hp, lgdtext, 'Location', 'NorthEast');
legend('boxoff');

result.zscore_bite = zscore_bite_all;
result.zscore_firstbite = zscore_firstbite_all;
result.t = t_aligned;
assignin('base', 'result', result);