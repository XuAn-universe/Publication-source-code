function AverageOrientation(app, Exp_Path, FrameRate, nmedian)
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
orientationxy_all_NI = cell(0);
orientationxychange_all_NI = cell(0);
frames_NI = 0;
n_I = 0;
orientationxy_all_I = cell(0);
orientationxychange_all_I = cell(0);
frames_I = 0;

laserstart_all = [];
laserstop_all = [];

orientationxybite_NI = [];
orientationxybite_I = [];
orientationxychangebite_NI = [];
orientationxychangebite_I = [];

orientationxyadj_NI = [];
orientationxyadj_I = [];

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

nadj_NI = zeros(1, ninhibition);
nadj_I = zeros(1, ninhibition);

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
            %             [~, ~, ~, ~, ~, x_top, y_top, z_top, ~, ~, laserstart, laserstop] = trajectory_postprocessing(21, Exp_Paths{i}, j,...
            %                 label_table, nmedian, FrameRate);
            %             [~, ~, ~, ~, ~, x_bottom, y_bottom, z_bottom, ~, ~, ~, ~] = trajectory_postprocessing(22, Exp_Paths{i}, j,...
            %                 label_table, nmedian, FrameRate);
            %             [~, ~, ~, ~, ~, x_center, y_center, z_center, ~, ~, ~, ~] = trajectory_postprocessing(33, Exp_Paths{i}, j,...
            %                 label_table, nmedian, FrameRate);
            %             x = [x_top(:, 1) x_bottom(:, 1) x_center(:, 1)];
            %             y = [y_top(:, 1) y_bottom(:, 1) y_center(:, 1)];
            %             z = [z_top(:, 1) z_bottom(:, 1) z_center(:, 1)];
            %             orientationxy = nan(size(x));
            %             orientationxy(:, 1) = atand((z(:, 1)-z(:, 2))./(sqrt((x(:, 1)-x(:, 2)).^2+(y(:, 1)-y(:, 2)).^2)));
            %             orientationxy(:, 2) = atand((z(:, 1)-z(:, 3))./(sqrt((x(:, 1)-x(:, 3)).^2+(y(:, 1)-y(:, 3)).^2)));
            %             orientationxy(:, 3) = atand((z(:, 2)-z(:, 3))./(sqrt((x(:, 2)-x(:, 3)).^2+(y(:, 2)-y(:, 3)).^2)));
            %             orientationxy = abs(orientationxy);
            %             orientationxy = mean(orientationxy, 2, 'omitnan');
            
            events_related2drop = [];
            adjustment_start = [];
            try
                temp = load([Exp_Paths{i} '\LabelledEvents' num2str(j) '.mat']);
                LabelledEvents = temp.LabelledEvents;
                MouthOpen = LabelledEvents.MouthOpen;
                TongueOut = LabelledEvents.TongueOut;
                PawLReachStart = LabelledEvents.PawLReachStart;
                PawRReachStart = LabelledEvents.PawRReachStart;
                events_related2drop = [MouthOpen; TongueOut; PawLReachStart; PawRReachStart];
                
                PawRAdjustmentStart = LabelledEvents.PawRAdjustmentStart;
                PawLAdjustmentStart = LabelledEvents.PawLAdjustmentStart;
                PawRAdjustmentEnd = LabelledEvents.PawRAdjustmentEnd;
                PawLAdjustmentEnd = LabelledEvents.PawLAdjustmentEnd;
                adjustment_start = get_adjustment_start(PawRAdjustmentStart, PawRAdjustmentEnd, PawLAdjustmentStart, PawLAdjustmentEnd);
            catch
                errordlg([Exp_Paths{i} '\LabelledEvents' num2str(j) '.mat is missing!'], 'Error');
            end
            
            if app.sDelay4sOnButton.Value
                if any(events_related2drop >= 4 & events_related2drop <= 8) % mouse drops pasta during inhibition
                    continue;
                end
            end

            [y_topPG1, y_topPG3, z_topPG1, z_topPG2, z_topPG3, x_top, ~, ~, ~, ~, laserstart, laserstop] = trajectory_postprocessing(21, Exp_Paths{i}, j,...
                label_table, nmedian, FrameRate);
            [y_bottomPG1, y_bottomPG3, z_bottomPG1, z_bottomPG2, z_bottomPG3, x_bottom, ~, ~, ~, ~, ~, ~] = trajectory_postprocessing(22, Exp_Paths{i}, j,...
                label_table, nmedian, FrameRate);
            [y_centerPG1, y_centerPG3, z_centerPG1, z_centerPG2, z_centerPG3, x_center, ~, ~, ~, ~, ~, ~] = trajectory_postprocessing(33, Exp_Paths{i}, j,...
                label_table, nmedian, FrameRate);
            x_all = [x_top(:, 1) x_center(:, 1) x_bottom(:, 1)];
            y_allPG1 = [y_topPG1(:, 1) y_centerPG1(:, 1) y_bottomPG1(:, 1)];
            y_allPG3 = [y_topPG3(:, 1) y_centerPG3(:, 1) y_bottomPG3(:, 1)];
            z_allPG1 = [z_topPG1(:, 1) z_centerPG1(:, 1) z_bottomPG1(:, 1)];
            z_allPG2 = [z_topPG2(:, 1) z_centerPG2(:, 1) z_bottomPG2(:, 1)];
            z_allPG3 = [z_topPG3(:, 1) z_centerPG3(:, 1) z_bottomPG3(:, 1)];
            avc = slope_estimation(x_all, z_allPG2);
%             avc = mean(avc, 2, 'omitnan');
            bvcPG1 = slope_estimation(y_allPG1, z_allPG1);
            bvcPG3 = slope_estimation(y_allPG3, z_allPG3);
            bvc = mean([bvcPG1 bvcPG3], 2, 'omitnan');
            orientationxy = acosd(1./sqrt(avc.^2+bvc.^2+1));
            orientationxy = 90-orientationxy;
            orientationxychange = [NaN; abs(diff(orientationxy))*FrameRate];

%             orientationxy_filtered = medfilt1(orientationxy, nmedian, 'omitnan', 'truncate');
            t = (1:size(orientationxy, 1))'/FrameRate;
            
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
                orientationxy_all_I{n_I} = orientationxy;
                orientationxychange_all_I{n_I} = orientationxychange;
                frames_I = max(frames_I, size(orientationxy, 1));
                
                if ~isempty(bite_timestamps)
                    for k = 1:ninhibition
                        for m = 1:numel(bite_timestamps)
                            timestamp = bite_timestamps(m);
                            if timestamp >= laser_timestamps(2*k-1) && timestamp <= laser_timestamps(2*k)
                                [~, id] = min(abs(t-timestamp));
                                if id-trange >= 1
                                    nbite_I(k) = nbite_I(k)+1;
                                    orientationxybite_I(:, nbite_I(k), k) = orientationxy(id-trange:id+trange);
                                    orientationxychangebite_I(:, nbite_I(k), k) = orientationxychange(id-trange:id+trange);
                                end
                            end
                        end
                    end
                end
                
                if ~isempty(adjustment_start)
                    for k = 1:ninhibition
                        for m = 1:numel(adjustment_start)
                            if 0
                                if (i == 1 && j == 2 && m == 8) || (i == 1 && j == 2 && m == 21)
                                    continue;
                                end
                            end
                                                
                            timestamp = adjustment_start(m);
                            if timestamp >= laser_timestamps(2*k-1) && timestamp <= laser_timestamps(2*k)
                                [~, id] = min(abs(t-timestamp));
                                if id-trange >= 1
                                    nadj_I(k) = nadj_I(k)+1;
                                    orientationxyadj_I(:, nadj_I(k), k) = orientationxy(id-trange:id+trange);
                                end
                            end
                        end
                    end
                end
            else
                n_NI = n_NI+1;
                orientationxy_all_NI{n_NI} = orientationxy;
                orientationxychange_all_NI{n_NI} = orientationxychange;
                frames_NI = max(frames_NI, size(orientationxy, 1));
                
                if ~isempty(bite_timestamps)
                    for k = 1:ninhibition
                        for m = 1:numel(bite_timestamps)
                            timestamp = bite_timestamps(m);
                            if timestamp >= tinhibition(2*k-1) && timestamp <= tinhibition(2*k)
                                [~, id] = min(abs(t-timestamp));
                                if id-trange >= 1
                                    nbite_NI(k) = nbite_NI(k)+1;
                                    orientationxybite_NI(:, nbite_NI(k), k) = orientationxy(id-trange:id+trange);
                                    orientationxychangebite_NI(:, nbite_NI(k), k) = orientationxychange(id-trange:id+trange);
                                end
                            end
                        end
                    end
                end
                
                if ~isempty(adjustment_start)
                    for k = 1:ninhibition
                        for m = 1:numel(adjustment_start)
                            if 0
                                if (i == 1 && j == 2 && m == 8) || (i == 1 && j == 2 && m == 21)
                                    continue;
                                end
                            end
                            
                            timestamp = adjustment_start(m);
                            if timestamp >= tinhibition(2*k-1) && timestamp <= tinhibition(2*k)
                                [~, id] = min(abs(t-timestamp));
                                if id-trange >= 1
                                    nadj_NI(k) = nadj_NI(k)+1;
                                    orientationxyadj_NI(:, nadj_NI(k), k) = orientationxy(id-trange:id+trange);
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end

orientationxy_matrix_NI = NaN(frames_NI, n_NI);
orientationxychange_matrix_NI = NaN(frames_NI, n_NI);

for i = 1:ninhibition
    orientationxylight_all_NI{i} = [];
    
    orientationxylight_mean_NI{i} = [];
    
    orientationxylight_var_NI{i} = [];
    
    orientationxychangelight_all_NI{i} = [];
    
    orientationxychangelight_mean_NI{i} = [];
end
for i = 1:n_NI
    frames_temp = size(orientationxy_all_NI{i}, 1);
    orientationxy_matrix_NI(1:frames_temp, i) = orientationxy_all_NI{i};
    orientationxychange_matrix_NI(1:frames_temp, i) = orientationxychange_all_NI{i};
    
    t = (1:frames_temp)'/FrameRate;
    for j = 1:ninhibition
        orientationxylight_all_NI{j} = [orientationxylight_all_NI{j}; orientationxy_all_NI{i}(t >= tinhibition(2*j-1) & t <= tinhibition(2*j))];
        
        orientationxylight_mean_NI{j} = [orientationxylight_mean_NI{j}; mean(orientationxy_all_NI{i}(t >= tinhibition(2*j-1) & t <= tinhibition(2*j)), 'omitnan')];
        
        orientationxylight_var_NI{j} = [orientationxylight_var_NI{j}; var(orientationxy_all_NI{i}(t >= tinhibition(2*j-1) & t <= tinhibition(2*j)), 'omitnan')];
        
        orientationxychangelight_all_NI{j} = [orientationxychangelight_all_NI{j}; orientationxychange_all_NI{i}(t >= tinhibition(2*j-1) & t <= tinhibition(2*j))];
        
        orientationxychangelight_mean_NI{j} = [orientationxychangelight_mean_NI{j}; mean(orientationxychange_all_NI{i}(t >= tinhibition(2*j-1) & t <= tinhibition(2*j)), 'omitnan')];
    end
end
t_NI = (1:frames_NI)'/FrameRate;
time_range = [app.timerangesEditField.Value app.toEditField.Value];

if ~app.NoLightButton.Value
    orientationxy_matrix_I = NaN(frames_I, n_I);
    orientationxychange_matrix_I = NaN(frames_I, n_I);
    
    for i = 1:ninhibition
        orientationxylight_all_I{i} = [];
        
        orientationxylight_mean_I{i} = [];
        
        orientationxylight_var_I{i} = [];
        
        orientationxychangelight_all_I{i} = [];
        
        orientationxychangelight_mean_I{i} = [];
    end
    for i = 1:n_I
        frames_temp = size(orientationxy_all_I{i}, 1);
        orientationxy_matrix_I(1:frames_temp, i) = orientationxy_all_I{i};
        orientationxychange_matrix_I(1:frames_temp, i) = orientationxychange_all_I{i};
        
        for j = 1:ninhibition
            orientationxylight_all_I{j} = [orientationxylight_all_I{j}; orientationxy_all_I{i}(laserstart_all(j, i):laserstop_all(j, i))];
            
            orientationxylight_mean_I{j} = [orientationxylight_mean_I{j}; mean(orientationxy_all_I{i}(laserstart_all(j, i):laserstop_all(j, i)), 'omitnan')];
            
            orientationxylight_var_I{j} = [orientationxylight_var_I{j}; var(orientationxy_all_I{i}(laserstart_all(j, i):laserstop_all(j, i)), 'omitnan')];
            
            orientationxychangelight_all_I{j} = [orientationxychangelight_all_I{j}; orientationxychange_all_I{i}(laserstart_all(j, i):laserstop_all(j, i))];
            
            orientationxychangelight_mean_I{j} = [orientationxychangelight_mean_I{j}; mean(orientationxychange_all_I{i}(laserstart_all(j, i):laserstop_all(j, i)), 'omitnan')];
        end
    end
    laserstart = mean(laserstart_all, 2)/FrameRate;
    laserstop = mean(laserstop_all, 2)/FrameRate;
    t_I = (1:frames_I)'/FrameRate;
    
    figure;
    plot_tj_multitrial(time_range, t_NI, orientationxy_matrix_NI, t_I, orientationxy_matrix_I,...
        laserstart, laserstop, ['Orientation (' char(176) ')'], 'REF XY Plane');
    figure;
    plot_tj_multitrial(time_range, t_NI, orientationxychange_matrix_NI, t_I, orientationxychange_matrix_I,...
        laserstart, laserstop, ['Orientation (' char(176) '/s)'], 'REF XY Plane');
else
    figure;
    plot_tj_multitrial(time_range, t_NI, orientationxy_matrix_NI, [], [],...
        laserstart, laserstop, ['Orientation (' char(176) ')'], 'REF XY Plane');
    figure;
    plot_tj_multitrial(time_range, t_NI, orientationxychange_matrix_NI, [], [],...
        laserstart, laserstop, ['Orientation (' char(176) '/s)'], 'REF XY Plane');
end

t = (-trange:1:trange)'/FrameRate;
time_range = [t(1) t(end)];
if ~app.NoLightButton.Value
    for i = 1:ninhibition
        % trajectories aligned to bite time
        figure;
        plot_tj_multitrial(time_range, t, orientationxybite_NI(:, 1:nbite_NI(i), i), t, orientationxybite_I(:, 1:nbite_I(i), i), 0, NaN, ['Orientation (' char(176) ')'], ['Bite, REF XY Plane (Inhibition #' num2str(i) ')']);
        figure;
        plot_tj_multitrial(time_range, t, orientationxychangebite_NI(:, 1:nbite_NI(i), i), t, orientationxychangebite_I(:, 1:nbite_I(i), i), 0, NaN, ['Orientation (' char(176) '/s)'], ['REF XY Plane (Inhibition #' num2str(i) ')']);
        figure;
        plot_tj_multitrial(time_range, t, orientationxyadj_NI(:, 1:nadj_NI(i), i), t, orientationxyadj_I(:, 1:nadj_I(i), i), 0, NaN, ['Orientation (' char(176) ')'], ['Adj, REF XY Plane (Inhibition #' num2str(i) ')']);
        compare2distributions(orientationxyadj_NI(trange+1, 1:nadj_NI(i), i), orientationxyadj_I(trange+1, 1:nadj_I(i), i), ['Orientation (' char(176) ')'], ['Adj, REF XY Plane (Inhibition #' num2str(i) ')']);
        
        % histograms
        if ~isempty(orientationxybite_I)
            histogram_multitrial(orientationxylight_all_NI{i}, orientationxybite_NI(trange+1, 1:nbite_NI(i), i), orientationxylight_all_I{i}, orientationxybite_I(trange+1, 1:nbite_I(i), i), ['Orientation (' char(176) ')'], ['REF XY Plane (Inhibition #' num2str(i) ')']);
            histogram_multitrial(orientationxychangelight_all_NI{i}, orientationxychangebite_NI(trange+1, 1:nbite_NI(i), i), orientationxychangelight_all_I{i}, orientationxychangebite_I(trange+1, 1:nbite_I(i), i), ['Orientation (' char(176) '/s)'], ['REF XY Plane (Inhibition #' num2str(i) ')']);
        else
            histogram_multitrial(orientationxylight_all_NI{i}, orientationxybite_NI(trange+1, 1:nbite_NI(i), i), orientationxylight_all_I{i}, [], ['Orientation (' char(176) ')'], ['REF XY Plane (Inhibition #' num2str(i) ')']);
            histogram_multitrial(orientationxychangelight_all_NI{i}, orientationxychangebite_NI(trange+1, 1:nbite_NI(i), i), orientationxychangelight_all_I{i}, [], ['Orientation (' char(176) '/s)'], ['REF XY Plane (Inhibition #' num2str(i) ')']);
        end
        
        histogram_multitrial(log(orientationxylight_var_NI{i}), [], log(orientationxylight_var_I{i}), [], ['Log Orientation Variance (' char(176) '^{2})'], ['REF XY Plane (Inhibition #' num2str(i) ')']);
        
        result.orientationxylight_all_NI(i) = orientationxylight_all_NI(i);
        
        result.orientationxylight_mean_NI(i) = mean(orientationxylight_mean_NI{i}, 'omitnan');
        
        result.orientationxylight_var_NI(i) = mean(orientationxylight_var_NI{i}, 'omitnan');
        
        result.orientationxybite_NI{i} = orientationxybite_NI(:, 1:nbite_NI(i), i);
        
        result.orientationxyadj_NI{i} = orientationxyadj_NI(:, 1:nadj_NI(i), i);
        
        result.orientationxychangelight_all_NI(i) = orientationxychangelight_all_NI(i);
        
        result.orientationxychangelight_mean_NI(i) = mean(orientationxychangelight_mean_NI{i}, 'omitnan');
        
        result.orientationxychangebite_NI{i} = orientationxychangebite_NI(:, 1:nbite_NI(i), i);
        
        result.orientationxylight_mean_I(i) = mean(orientationxylight_mean_I{i}, 'omitnan');
        
        result.orientationxylight_var_I(i) = mean(orientationxylight_var_I{i}, 'omitnan');
        
        result.orientationxylight_all_I(i) = orientationxylight_all_I(i);
        
        result.orientationxybite_I{i} = orientationxybite_I(:, 1:nbite_I(i), i);
        
        result.orientationxyadj_I{i} = orientationxyadj_I(:, 1:nadj_I(i), i);
        
        result.orientationxychangelight_all_I(i) = orientationxychangelight_all_I(i);
        
        result.orientationxychangelight_mean_I(i) = mean(orientationxychangelight_mean_I{i}, 'omitnan');
        
        result.orientationxychangebite_I{i} = orientationxychangebite_I(:, 1:nbite_I(i), i);
    end
else
    for i = 1:ninhibition
        % trajectories aligned to bite time
        figure;
        plot_tj_multitrial(time_range, t, orientationxybite_NI(:, 1:nbite_NI(i), i), [], [], 0, NaN, ['Orientation (' char(176) ')'], ['Bite, REF XY Plane (Inhibition #' num2str(i) ')']);
        figure;
        plot_tj_multitrial(time_range, t, orientationxychangebite_NI(:, 1:nbite_NI(i), i), [], [], 0, NaN, ['Orientation (' char(176) '/s)'], ['REF XY Plane (Inhibition #' num2str(i) ')']);
        figure;
        plot_tj_multitrial(time_range, t, orientationxyadj_NI(:, 1:nadj_NI(i), i), [], [], 0, NaN, ['Orientation (' char(176) ')'], ['Adj, REF XY Plane (Inhibition #' num2str(i) ')']);
        
        % histograms
        histogram_multitrial(orientationxylight_all_NI{i}, orientationxybite_NI(trange+1, 1:nbite_NI(i), i), [], [], ['Orientation (' char(176) ')'], ['REF XY Plane (Inhibition #' num2str(i) ')']);
        histogram_multitrial(orientationxychangelight_all_NI{i}, orientationxychangebite_NI(trange+1, 1:nbite_NI(i), i), [], [], ['Orientation (' char(176) '/s)'], ['REF XY Plane (Inhibition #' num2str(i) ')']);
        
        histogram_multitrial(log(orientationxylight_var_NI{i}), [], [], [], ['Log Orientation Variance (' char(176) '^{2})'], ['REF XY Plane (Inhibition #' num2str(i) ')']);
        
        result.orientationxylight_all_NI(i) = orientationxylight_all_NI(i);
        
        result.orientationxylight_mean_NI(i) = mean(orientationxylight_mean_NI{i}, 'omitnan');
        
        result.orientationxylight_var_NI(i) = mean(orientationxylight_var_NI{i}, 'omitnan');
        
        result.orientationxybite_NI{i} = orientationxybite_NI(:, 1:nbite_NI(i), i);
        
        result.orientationxyadj_NI{i} = orientationxyadj_NI(:, 1:nadj_NI(i), i);
        
        result.orientationxychangelight_all_NI(i) = orientationxychangelight_all_NI(i);
        
        result.orientationxychangelight_mean_NI(i) = mean(orientationxychangelight_mean_NI{i}, 'omitnan');
        
        result.orientationxychangebite_NI{i} = orientationxychangebite_NI(:, 1:nbite_NI(i), i);
    end
end

[file, path] = uiputfile('Pasta Orientation.mat', 'Save result');
if file ~= 0
    save([path file], 'result');
end