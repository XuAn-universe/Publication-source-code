function TrajectorytoEvents(app, Exp_Path, FrameRate, nmedian)
persistent path
value = app.TrialsListBox.Value;
value = sort(value);
trials = numel(value);

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
trange = app.bitesEditField.Value; % time range for event alignment
trange = trange*FrameRate;

zbite_NI = [];
zbite_NI_1p = cell(1, ninhibition);
zspeedbite_NI = [];
xyzspeedbite_NI = [];
zaccelerationbite_NI = [];
xyzaccelerationbite_NI = [];
zsit_NI = [];
zsit_NI_1p = [];
zpawlreach_NI = [];
zpawlreach_NI_1p = [];
zpawrreach_NI = [];
zpawrreach_NI_1p = [];
dzwithdraw_NI = [];
dzspeedwithdraw_NI = [];

zbite_I = [];
zbite_I_1p = cell(1, ninhibition);
zspeedbite_I = [];
xyzspeedbite_I = [];
zaccelerationbite_I = [];
xyzaccelerationbite_I = [];
zsit_I = [];
zsit_I_1p = [];
zpawlreach_I = [];
zpawlreach_I_1p = [];
zpawrreach_I = [];
zpawrreach_I_1p = [];
dzwithdraw_I = [];
dzspeedwithdraw_I = [];

nbite_NI = zeros(1, ninhibition);
nbite_I = zeros(1, ninhibition);
nsit_NI = 0;
nsit_I = 0;
npawlreach_NI = 0;
npawlreach_I = 0;
npawrreach_NI = 0;
npawrreach_I = 0;
nwithdraw_NI = zeros(1, ninhibition);
nwithdraw_I = zeros(1, ninhibition);

% process trajectories
if app.TrackingDataCheckBox.Value
    label_table = table2array(app.UITable.Data);
else
    helpdlg('''Tracking Data'' needs to be checked to proceed');
    return;
end
pointID = app.DropDown.Value;

% get bite events
Bite_events = [];
try
    audiolocation = Exp_Path(1:end-7);
    temp = load([audiolocation '\Detected_Bite_Events.mat']);
    Bite_events = temp.Audio_analysis;
end

for i = 1:trials
    PawLReachEnd = [];
    PawRReachEnd = [];
    SitEnd = [];
    Withdraw = [];
    try
        temp = load([Exp_Path '\LabelledEvents' num2str(value(i)) '.mat']);
        LabelledEvents = temp.LabelledEvents;
        PawLReachEnd = LabelledEvents.PawLReachEnd;
        PawRReachEnd = LabelledEvents.PawRReachEnd;
        SitEnd = LabelledEvents.SitEnd;
        Withdraw = LabelledEvents.BiteBoutStart;
    end
    
    bite_timestamps = [];
    if ~isempty(Bite_events)
        bite_timestamps = Bite_events(value(i)).time_bites;
        if ~isempty(bite_timestamps)
            bite_amplitudes = Bite_events(value(i)).amplitude_bites/(max(Bite_events(value(i)).amplitude_bites));
        end
    end
    
    [~, ~, ~, ~, ~, x, y, z, speed, acceleration, laserstart, laserstop] = trajectory_postprocessing(pointID, Exp_Path, value(i), label_table, nmedian, FrameRate);
    [~, ~, ~, ~, ~, x_nose, y_nose, z_nose, speed_nose, acceleration_nose, ~, ~] = trajectory_postprocessing(3, Exp_Path, value(i), label_table, nmedian, FrameRate);
    dz = z_nose-z;
    dz = dz(:, 1);
    dzspeed = diff(dz)*FrameRate;
    dzspeed = [NaN; dzspeed];
    t = (1:size(x, 1))'/FrameRate;
    
    if app.NoLightButton.Value
        laserstart = [];
        laserstop = [];
    else
        laser_timestamps = Bite_events(value(i)).laser_timestamps;
    end
    if ~isempty(laserstart)
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
                        end
                    end
                end
            end
        end
        
        if ~isempty(Withdraw)
            for k = 1:ninhibition
                for m = 1:numel(Withdraw)
                    timestamp = Withdraw(m);
                    if timestamp >= laser_timestamps(2*k-1) && timestamp <= laser_timestamps(2*k)
                        [~, id] = min(abs(t-timestamp));
                        if id-trange >= 1
                            nwithdraw_I(k) = nwithdraw_I(k)+1;
                            dzwithdraw_I(:, nwithdraw_I(k), k) = dz(id-trange:id+trange, 1);
                            dzspeedwithdraw_I(:, nwithdraw_I(k), k) = abs(dzspeed(id-trange:id+trange, 1));
                        end
                    end
                end
            end
        end
    else
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
                        end
                    end
                end
            end
        end
        
        if ~isempty(Withdraw)
            for k = 1:ninhibition
                for m = 1:numel(Withdraw)
                    timestamp = Withdraw(m);
                    if timestamp >= tinhibition(2*k-1) && timestamp <= tinhibition(2*k)
                        [~, id] = min(abs(t-timestamp));
                        if id-trange >= 1
                            nwithdraw_NI(k) = nwithdraw_NI(k)+1;
                            dzwithdraw_NI(:, nwithdraw_NI(k), k) = dz(id-trange:id+trange, 1);
                            dzspeedwithdraw_NI(:, nwithdraw_NI(k), k) = abs(dzspeed(id-trange:id+trange, 1));
                        end
                    end
                end
            end
        end
        
        if app.NoLightButton.Value
            if pointID == 3 && ~isempty(SitEnd)
                timestamp = SitEnd(1);
                if timestamp >= tinhibition(1) && timestamp <= tinhibition(2)
                    [~, id] = min(abs(t-timestamp));
                    if id-trange >= 1
                        nsit_NI = nsit_NI+1;
                        zsit_NI(:, nsit_NI) = z(id-trange:id+trange, 1);
                    end
                end
            end
            if pointID == 10 && ~isempty(PawLReachEnd)
                timestamp = PawLReachEnd(1);
                if timestamp >= tinhibition(1) && timestamp <= tinhibition(2)
                    [~, id] = min(abs(t-timestamp));
                    if id-trange >= 1
                        npawlreach_NI = npawlreach_NI+1;
                        zpawlreach_NI(:, npawlreach_NI) = z(id-trange:id+trange, 1);
                    end
                end
            end
            if pointID == 16 && ~isempty(PawRReachEnd)
                timestamp = PawRReachEnd(1);
                if timestamp >= tinhibition(1) && timestamp <= tinhibition(2)
                    [~, id] = min(abs(t-timestamp));
                    if id-trange >= 1
                        npawrreach_NI = npawrreach_NI+1;
                        zpawrreach_NI(:, npawrreach_NI) = z(id-trange:id+trange, 1);
                    end
                end
            end
        end
    end
end

t = (-trange:1:trange)'/FrameRate;
time_range = [t(1) t(end)];
if ~app.NoLightButton.Value
    for i = 1:ninhibition
        % trajectories aligned to withdraw time
        try
            figure;
            plot_tj_multitrial(time_range, t, dzwithdraw_NI(:, 1:nwithdraw_NI(i), i), t, dzwithdraw_I(:, 1:nwithdraw_I(i), i), 0, NaN, 'dZ (mm)', ['Inhibition #' num2str(i)]);
            figure;
            plot_tj_multitrial(time_range, t, dzspeedwithdraw_NI(:, 1:nwithdraw_NI(i), i), t, dzspeedwithdraw_I(:, 1:nwithdraw_I(i), i), 0, NaN, 'Speed (mm/s)', ['dZ (Inhibition #' num2str(i) ')']);
        end
        
        result.dzwithdraw_NI{i} = dzwithdraw_NI(:, 1:nwithdraw_NI(i), i);
        result.dzspeedwithdraw_NI{i} = dzspeedwithdraw_NI(:, 1:nwithdraw_NI(i), i);
        
        result.dzwithdraw_I{i} = dzwithdraw_I(:, 1:nwithdraw_I(i), i);
        result.dzspeedwithdraw_I{i} = dzspeedwithdraw_I(:, 1:nwithdraw_I(i), i);
    end
else
    if ~isempty(dzwithdraw_NI)
        for i = 1:ninhibition
            % trajectories aligned to bite time
            figure;
            plot_tj_multitrial(time_range, t, dzwithdraw_NI(:, 1:nwithdraw_NI(i), i), [], [], 0, NaN, 'dZ (mm)', ['Inhibition #' num2str(i)]);
            figure;
            plot_tj_multitrial(time_range, t, dzspeedwithdraw_NI(:, 1:nwithdraw_NI(i), i), [], [], 0, NaN, 'Speed (mm/s)', ['dZ (Inhibition #' num2str(i) ')']);
            
            result.dzwithdraw_NI{i} = dzwithdraw_NI(:, 1:nwithdraw_NI(i), i);
            result.dzspeedwithdraw_NI{i} = dzspeedwithdraw_NI(:, 1:nwithdraw_NI(i), i);
        end
    end
end

if ~app.NoLightButton.Value
    try
        for i = 1:ninhibition
            zbite_NI_1p{i} = zbite_NI(trange+1, 1:nbite_NI(i), i);
            if nbite_I(i) ~= 0
                zbite_I_1p{i} = zbite_I(trange+1, 1:nbite_I(i), i);
            end
        end
        if ~isempty(zsit_NI) && ~isempty(zsit_I)
            zsit_NI_1p = zsit_NI(trange+1, :);
            zsit_I_1p = zsit_I(trange+1, :);
        end
        if ~isempty(zpawlreach_NI) && ~isempty(zpawlreach_I)
            zpawlreach_NI_1p = zpawlreach_NI(trange+1, :);
            zpawlreach_I_1p = zpawlreach_I(trange+1, :);
        end
        if ~isempty(zpawrreach_NI) && ~isempty(zpawrreach_I)
            zpawrreach_NI_1p = zpawrreach_NI(trange+1, :);
            zpawrreach_I_1p = zpawrreach_I(trange+1, :);
        end
    end
else
    for i = 1:ninhibition
        zbite_NI_1p{i} = zbite_NI(trange+1, 1:nbite_NI(i), i);
    end
    result.zbite_NI_1p = zbite_NI_1p;
    if ~isempty(zsit_NI)
        zsit_NI_1p = zsit_NI(trange+1, :);
        result.zsit_NI_1p = zsit_NI_1p;
    end
    if ~isempty(zpawlreach_NI)
        zpawlreach_NI_1p = zpawlreach_NI(trange+1, :);
        result.zpawlreach_NI_1p = zpawlreach_NI_1p;
    end
    if ~isempty(zpawrreach_NI)
        zpawrreach_NI_1p = zpawrreach_NI(trange+1, :);
        result.zpawrreach_NI_1p = zpawrreach_NI_1p;
    end
end
assignin('base', 'result', result);

curpwd = pwd;
try
    cd(path);
end
[file, path] = uiputfile('Zposition4events.mat', 'Save result');
if file ~= 0
    save([path file], 'result');
end
cd(curpwd);