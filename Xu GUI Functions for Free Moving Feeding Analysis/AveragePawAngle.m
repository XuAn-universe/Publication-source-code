function AveragePawAngle(app, Exp_Path, FrameRate, nmedian, pawsegment, xyzspeed_threshold)
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

n_NI = 0;
pawangle_all_NI = cell(0);
frames_NI = 0;

n_I = 0;
pawangle_all_I = cell(0);
frames_I = 0;
laserstart_all = [];
laserstop_all = [];

pawanglebite_NI = [];
pawanglebite_I = [];

if app.sOnTwiceButton.Value
    ninhibition = 2;
    tinhibition = [4 13; 8 17];
elseif app.NoLightButton.Value || app.WholeTrialOnButton.Value
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
            [~, ~, ~, ~, ~, x_PL, y_PL, z_PL, ~, ~, laserstart, laserstop] = trajectory_postprocessing(10, Exp_Paths{i}, j,...
                label_table, nmedian, FrameRate);
            [~, ~, ~, ~, ~, x_PR, y_PR, z_PR, ~, ~, ~, ~] = trajectory_postprocessing(16, Exp_Paths{i}, j,...
                label_table, nmedian, FrameRate);
            vPL = [x_PL(:, 1) y_PL(:, 1) z_PL(:, 1)];
            vPL = circshift(vPL, -pawsegment, 1)-vPL;
            vPL(end-pawsegment+1:end, :) = nan;
            vPR = [x_PR(:, 1) y_PR(:, 1) z_PR(:, 1)];
            vPR = circshift(vPR, -pawsegment, 1)-vPR;
            vPR(end-pawsegment+1:end, :) = nan;
            vangle = atan2d(vecnorm(cross(vPL, vPR, 2), 2, 2), dot(vPL, vPR, 2));
%             vangle_filtered = medfilt1(vangle, nmedian, 'omitnan', 'truncate');

            t = (1:size(vangle, 1))'/FrameRate;
            
            bite_timestamps = Bite_events(j).time_bites;
            if app.NoLightButton.Value
                laserstart = [];
                laserstop = [];
            else
                laser_timestamps = Bite_events(j).laser_timestamps;
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
                pawangle_all_I{n_I} = vangle;
                frames_I = max(frames_I, size(vangle, 1));
                
                if ~isempty(bite_timestamps)
                    for k = 1:ninhibition
                        for m = 1:numel(bite_timestamps)
                            timestamp = bite_timestamps(m);
                            if timestamp >= laser_timestamps(2*k-1) && timestamp <= laser_timestamps(2*k)
                                [~, id] = min(abs(t-timestamp));
                                if id-trange >= 1
                                    nbite_I(k) = nbite_I(k)+1;
                                    pawanglebite_I(:, nbite_I(k), k) = vangle(id-trange:id+trange);
                                end
                            end
                        end
                    end
                end
            else
                n_NI = n_NI+1;
                pawangle_all_NI{n_NI} = vangle;
                frames_NI = max(frames_NI, size(vangle, 1));
                
                if ~isempty(bite_timestamps)
                    for k = 1:ninhibition
                        for m = 1:numel(bite_timestamps)
                            timestamp = bite_timestamps(m);
                            if timestamp >= tinhibition(2*k-1) && timestamp <= tinhibition(2*k)
                                [~, id] = min(abs(t-timestamp));
                                if id-trange >= 1
                                    nbite_NI(k) = nbite_NI(k)+1;
                                    pawanglebite_NI(:, nbite_NI(k), k) = vangle(id-trange:id+trange);
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end

pawangle_matrix_NI = NaN(frames_NI, n_NI);

for i = 1:ninhibition
    pawanglelight_all_NI{i} = [];
    
    pawanglelight_mean_NI{i} = [];
    
    pawanglelight_var_NI{i} = [];
end
for i = 1:n_NI
    frames_temp = size(pawangle_all_NI{i}, 1);
    pawangle_matrix_NI(1:frames_temp, i) = pawangle_all_NI{i};
    
    t = (1:frames_temp)'/FrameRate;
    for j = 1:ninhibition
        pawanglelight_all_NI{j} = [pawanglelight_all_NI{j}; pawangle_all_NI{i}(t >= tinhibition(2*j-1) & t <= tinhibition(2*j))];
        
        pawanglelight_mean_NI{j} = [pawanglelight_mean_NI{j}; mean(pawangle_all_NI{i}(t >= tinhibition(2*j-1) & t <= tinhibition(2*j)), 'omitnan')];
        
        pawanglelight_var_NI{j} = [pawanglelight_var_NI{j}; var(pawangle_all_NI{i}(t >= tinhibition(2*j-1) & t <= tinhibition(2*j)), 'omitnan')];
    end
end
t_NI = (1:frames_NI)'/FrameRate;
time_range = [app.timerangesEditField.Value app.toEditField.Value];

if ~app.NoLightButton.Value
    pawangle_matrix_I = NaN(frames_I, n_I);
    
    for i = 1:ninhibition
        pawanglelight_all_I{i} = [];
        
        pawanglelight_mean_I{i} = [];
        
        pawanglelight_var_I{i} = [];
    end
    for i = 1:n_I
        frames_temp = size(pawangle_all_I{i}, 1);
        pawangle_matrix_I(1:frames_temp, i) = pawangle_all_I{i};
        
        for j = 1:ninhibition
            pawanglelight_all_I{j} = [pawanglelight_all_I{j}; pawangle_all_I{i}(laserstart_all(j, i):laserstop_all(j, i))];
            
            pawanglelight_mean_I{j} = [pawanglelight_mean_I{j}; mean(pawangle_all_I{i}(laserstart_all(j, i):laserstop_all(j, i)), 'omitnan')];
            
            pawanglelight_var_I{j} = [pawanglelight_var_I{j}; var(pawangle_all_I{i}(laserstart_all(j, i):laserstop_all(j, i)), 'omitnan')];
        end
    end
    laserstart = mean(laserstart_all, 2)/FrameRate;
    laserstop = mean(laserstop_all, 2)/FrameRate;
    t_I = (1:frames_I)'/FrameRate;
    
    figure;
    plot_tj_multitrial(time_range, t_NI, pawangle_matrix_NI, t_I, pawangle_matrix_I,...
        laserstart, laserstop, ['Angle (' char(176) ')'], 'Angle of Paw Vectors');
else
    figure;
    plot_tj_multitrial(time_range, t_NI, pawangle_matrix_NI, [], [],...
        laserstart, laserstop, ['Angle (' char(176) ')'], 'Angle of Paw Vectors');
end

t = (-trange:1:trange)'/FrameRate;
time_range = [t(1) t(end)];
if ~app.NoLightButton.Value
    for i = 1:ninhibition
        % trajectories aligned to bite time
        figure;
        plot_tj_multitrial(time_range, t, pawanglebite_NI(:, 1:nbite_NI(i), i), t, pawanglebite_I(:, 1:nbite_I(i), i), 0, NaN, ['Angle (' char(176) ')'], ['Angle of Paw Vectors (Inhibition #' num2str(i) ')']);
        
        % histograms
        histogram_multitrial(pawanglelight_all_NI{i}, pawanglebite_NI(trange+1, :, i), pawanglelight_all_I{i}, pawanglebite_I(trange+1, :, i), ['Angle (' char(176) ')'], ['Angle of Paw Vectors (Inhibition #' num2str(i) ')']);
        
        histogram_multitrial(log(pawanglelight_var_NI{i}), [], log(pawanglelight_var_I{i}), [], ['Log Angle Variance (' char(176) '^{2})'], ['Angle of Paw Vectors (Inhibition #' num2str(i) ')']);
        
        result.pawanglelight_mean_NI(i) = mean(pawanglelight_mean_NI{i}, 'omitnan');
        
        result.pawanglelight_var_NI(i) = mean(pawanglelight_var_NI{i}, 'omitnan');
        
        result.pawanglelight_mean_I(i) = mean(pawanglelight_mean_I{i}, 'omitnan');
        
        result.pawanglelight_var_I(i) = mean(pawanglelight_var_I{i}, 'omitnan');
    end
else
    for i = 1:ninhibition
        % trajectories aligned to bite time
        figure;
        plot_tj_multitrial(time_range, t, pawanglebite_NI(:, 1:nbite_NI(i), i), [], [], 0, NaN, ['Angle (' char(176) ')'], ['Angle of Paw Vectors (Inhibition #' num2str(i) ')']);
        
        % histograms
        histogram_multitrial(pawanglelight_all_NI{i}, pawanglebite_NI(trange+1, :, i), [], [], ['Angle (' char(176) ')'], ['Angle of Paw Vectors (Inhibition #' num2str(i) ')']);
        
        histogram_multitrial(log(pawanglelight_var_NI{i}), [], [], [], ['Log Angle Variance (' char(176) '^{2})'], ['Angle of Paw Vectors (Inhibition #' num2str(i) ')']);
        
        result.pawanglelight_mean_NI(i) = mean(pawanglelight_mean_NI{i}, 'omitnan');
        
        result.pawanglelight_var_NI(i) = mean(pawanglelight_var_NI{i}, 'omitnan');
    end
end

[file, path] = uiputfile('AveragePawAngle.mat', 'Save result');
if file ~= 0
    save([path file], 'result');
end