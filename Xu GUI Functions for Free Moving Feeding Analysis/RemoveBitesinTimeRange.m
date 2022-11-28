function RemoveBitesinTimeRange(app, Exp_Path, duration)
trial = app.TrialsListBox.Value;
if numel(trial) ~= 1
    helpdlg('You can only choose one trial');
    return;
end

% get bite events
audiolocation = Exp_Path(1:end-7);
try
    temp = load([audiolocation '\Detected_Bite_Events.mat']);
    Audio_analysis = temp.Audio_analysis;
    bite_timestamps = Audio_analysis(trial).time_bites;
    bite_amplitudes = Audio_analysis(trial).amplitude_bites;
    bite_amplitudes(bite_timestamps >= duration(1) & bite_timestamps <= duration(2)) = [];
    bite_timestamps(bite_timestamps >= duration(1) & bite_timestamps <= duration(2)) = [];
    Audio_analysis(trial).time_bites = bite_timestamps;
    Audio_analysis(trial).amplitude_bites = bite_amplitudes;
    save([audiolocation '\Detected_Bite_Events.mat'], 'Audio_analysis');
    msgbox('Done !');
catch
    errordlg('The bite events have not been detected yet!');
    return;
end