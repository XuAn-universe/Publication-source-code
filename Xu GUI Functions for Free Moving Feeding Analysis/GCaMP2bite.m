function GCaMP2bite(app, Exp_Path, fpdata)
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

if ~isempty(fpdata)
    nchannel = numel(fpdata);
    for i = 1:nchannel
        fpdata_zsignal(:, i) = fpdata{i}(:, 2);
        lgdtext{i} = ['Channel ' num2str(i)];
    end
    fpdata_t = fpdata{1}(:, 1);
    
    pre_bite = app.presEditField.Value;
    post_bite = app.postsEditField.Value;
    [t_bite, zscore_bite, zscore_bite_random] = AlignSignal2Event(fpdata_t, fpdata_zsignal, bite_timestamps, pre_bite, post_bite);
    for i = 1:nchannel
        figure;
        plot_gcamp_align2event(t_bite, zscore_bite_random(:, :, i), zscore_bite(:, :, i), ['aligned to bite: Channel ' num2str(i)]);
    end
end