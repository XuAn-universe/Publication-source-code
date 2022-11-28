function AverageDistance(app, Exp_Path, FrameRate, nmedian, zspeed_threshold)
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
pointID1 = app.DropDown_2.Value;
pointID2 = app.DropDown_3.Value;

n_NI = 0;
distance_all_NI = cell(0);
orientationxy_all_NI = cell(0);
dz1_all_NI = cell(0);
dz2_all_NI = cell(0);
z1_all_NI = cell(0);
z2_all_NI = cell(0);
frames_NI = 0;

n_I = 0;
distance_all_I = cell(0);
orientationxy_all_I = cell(0);
dz1_all_I = cell(0);
dz2_all_I = cell(0);
z1_all_I = cell(0);
z2_all_I = cell(0);
frames_I = 0;
laserstart_all = [];
laserstop_all = [];

dzbite_NI = [];
dxyzbite_NI = [];
orientationxybite_NI = [];

dzbite_I = [];
dxyzbite_I = [];
orientationxybite_I = [];

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
            
            [~, ~, ~, ~, ~, x1, y1, z1, ~, ~, laserstart, laserstop] = trajectory_postprocessing(pointID1, Exp_Paths{i}, j,...
                label_table, nmedian, FrameRate);
            [~, ~, ~, ~, ~, x2, y2, z2, ~, ~, ~, ~] = trajectory_postprocessing(pointID2, Exp_Paths{i}, j, label_table, nmedian, FrameRate);
            dx = x1(:, 1)-x2(:, 1);
%             dx_filtered = medfilt1(dx, nmedian, 'omitnan', 'truncate');
            dy = y1(:, 1)-y2(:, 1);
%             dy_filtered = medfilt1(dy, nmedian, 'omitnan', 'truncate');
            dz = z1(:, 1)-z2(:, 1);
%             dz_filtered = medfilt1(dz, nmedian, 'omitnan', 'truncate');
            dxy = sqrt(dx.^2+dy.^2);
%             dxy_filtered = medfilt1(dxy, nmedian, 'omitnan', 'truncate');
            dxz = sqrt(dx.^2+dz.^2);
%             dxz_filtered = medfilt1(dxz, nmedian, 'omitnan', 'truncate');
            dyz = sqrt(dy.^2+dz.^2);
%             dyz_filtered = medfilt1(dyz, nmedian, 'omitnan', 'truncate');
            dxyz = sqrt(dx.^2+dy.^2+dz.^2);
%             dxyz_filtered = medfilt1(dxyz, nmedian, 'omitnan', 'truncate');
            orientationxy = atand(dz./(sqrt(dx.^2+dy.^2)));
%             orientationxy_filtered = medfilt1(orientationxy, nmedian, 'omitnan', 'truncate');
            dz1 = [NaN; diff(z1(:, 1))];
%             dz1_filtered = medfilt1(dz1, nmedian, 'omitnan', 'truncate');
%             dz1_filtered(1) = NaN;
            dz2 = [NaN; diff(z2(:, 1))];
%             dz2_filtered = medfilt1(dz2, nmedian, 'omitnan', 'truncate');
%             dz2_filtered(1) = NaN;
            t = (1:size(dx, 1))'/FrameRate;
            
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
                
                distance_all_I{n_I}(:, 1) = dx;
                distance_all_I{n_I}(:, 2) = dy;
                distance_all_I{n_I}(:, 3) = dz;
                distance_all_I{n_I}(:, 4) = dxy;
                distance_all_I{n_I}(:, 5) = dxz;
                distance_all_I{n_I}(:, 6) = dyz;
                distance_all_I{n_I}(:, 7) = dxyz;
                orientationxy_all_I{n_I} = orientationxy;
                dz1_all_I{n_I} = dz1;
                dz2_all_I{n_I} = dz2;
                z1_all_I{n_I} = z1(:, 1);
                z2_all_I{n_I} = z2(:, 1);
                frames_I = max(frames_I, size(dx, 1));
                
                if ~isempty(bite_timestamps)
                    for k = 1:ninhibition
                        for m = 1:numel(bite_timestamps)
                            timestamp = bite_timestamps(m);
                            if timestamp >= laser_timestamps(2*k-1) && timestamp <= laser_timestamps(2*k)
                                [~, id] = min(abs(t-timestamp));
                                if id-trange >= 1
                                    nbite_I(k) = nbite_I(k)+1;
                                    dzbite_I(:, nbite_I(k), k) = dz(id-trange:id+trange);
                                    dxyzbite_I(:, nbite_I(k), k) = dxyz(id-trange:id+trange);
                                    orientationxybite_I(:, nbite_I(k), k) = orientationxy(id-trange:id+trange);
                                end
                            end
                        end
                    end
                end
            else
                n_NI = n_NI+1;
                distance_all_NI{n_NI}(:, 1) = dx;
                distance_all_NI{n_NI}(:, 2) = dy;
                distance_all_NI{n_NI}(:, 3) = dz;
                distance_all_NI{n_NI}(:, 4) = dxy;
                distance_all_NI{n_NI}(:, 5) = dxz;
                distance_all_NI{n_NI}(:, 6) = dyz;
                distance_all_NI{n_NI}(:, 7) = dxyz;
                orientationxy_all_NI{n_NI} = orientationxy;
                dz1_all_NI{n_NI} = dz1;
                dz2_all_NI{n_NI} = dz2;
                z1_all_NI{n_NI} = z1(:, 1);
                z2_all_NI{n_NI} = z2(:, 1);
                frames_NI = max(frames_NI, size(dx, 1));
                
                if ~isempty(bite_timestamps)
                    for k = 1:ninhibition
                        for m = 1:numel(bite_timestamps)
                            timestamp = bite_timestamps(m);
                            if timestamp >= tinhibition(2*k-1) && timestamp <= tinhibition(2*k)
                                [~, id] = min(abs(t-timestamp));
                                if id-trange >= 1
                                    nbite_NI(k) = nbite_NI(k)+1;
                                    dzbite_NI(:, nbite_NI(k), k) = dz(id-trange:id+trange);
                                    dxyzbite_NI(:, nbite_NI(k), k) = dxyz(id-trange:id+trange);
                                    orientationxybite_NI(:, nbite_NI(k), k) = orientationxy(id-trange:id+trange);
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end

distance_matrix_NI = NaN(frames_NI, 7, n_NI);
orientationxy_matrix_NI = NaN(frames_NI, n_NI);

for i = 1:ninhibition
    dzlight_all_NI{i} = [];
    dxyzlight_all_NI{i} = [];
    orientationxylight_all_NI{i} = [];
    
    dzlight_mean_NI{i} = [];
    dxyzlight_mean_NI{i} = [];
    orientationxylight_mean_NI{i} = [];
    
    dzlight_var_NI{i} = [];
    dxyzlight_var_NI{i} = [];
    orientationxylight_var_NI{i} = [];
    
    dz1_NI{i} = [];
    dz2_NI{i} = [];
    z1_NI{i} = [];
    z2_NI{i} = [];
end
for i = 1:n_NI
    frames_temp = size(distance_all_NI{i}, 1);
    distance_matrix_NI(1:frames_temp, :, i) = distance_all_NI{i};
    orientationxy_matrix_NI(1:frames_temp, i) = orientationxy_all_NI{i};
    
    t = (1:frames_temp)'/FrameRate;
    for j = 1:ninhibition
        dzlight_all_NI{j} = [dzlight_all_NI{j}; distance_all_NI{i}(t >= tinhibition(2*j-1) & t <= tinhibition(2*j), 3)];
        dxyzlight_all_NI{j} = [dxyzlight_all_NI{j}; distance_all_NI{i}(t >= tinhibition(2*j-1) & t <= tinhibition(2*j), 7)];
        orientationxylight_all_NI{j} = [orientationxylight_all_NI{j}; orientationxy_all_NI{i}(t >= tinhibition(2*j-1) & t <= tinhibition(2*j))];
        
        dzlight_mean_NI{j} = [dzlight_mean_NI{j}; mean(distance_all_NI{i}(t >= tinhibition(2*j-1) & t <= tinhibition(2*j), 3), 'omitnan')];
        dxyzlight_mean_NI{j} = [dxyzlight_mean_NI{j}; mean(distance_all_NI{i}(t >= tinhibition(2*j-1) & t <= tinhibition(2*j), 7), 'omitnan')];
        orientationxylight_mean_NI{j} = [orientationxylight_mean_NI{j}; mean(orientationxy_all_NI{i}(t >= tinhibition(2*j-1) & t <= tinhibition(2*j)), 'omitnan')];
        
        dzlight_var_NI{j} = [dzlight_var_NI{j}; var(distance_all_NI{i}(t >= tinhibition(2*j-1) & t <= tinhibition(2*j), 3), 'omitnan')];
        dxyzlight_var_NI{j} = [dxyzlight_var_NI{j}; var(distance_all_NI{i}(t >= tinhibition(2*j-1) & t <= tinhibition(2*j), 7), 'omitnan')];
        orientationxylight_var_NI{j} = [orientationxylight_var_NI{j}; var(orientationxy_all_NI{i}(t >= tinhibition(2*j-1) & t <= tinhibition(2*j)), 'omitnan')];
        
        dz1_NI{j} = [dz1_NI{j}; dz1_all_NI{i}(t >= tinhibition(2*j-1) & t <= tinhibition(2*j))];
        dz2_NI{j} = [dz2_NI{j}; dz2_all_NI{i}(t >= tinhibition(2*j-1) & t <= tinhibition(2*j))];
        z1_NI{j} = [z1_NI{j}; z1_all_NI{i}(t >= tinhibition(2*j-1) & t <= tinhibition(2*j))];
        z2_NI{j} = [z2_NI{j}; z2_all_NI{i}(t >= tinhibition(2*j-1) & t <= tinhibition(2*j))];
    end
end
t_NI = (1:frames_NI)'/FrameRate;
time_range = [app.timerangesEditField.Value app.toEditField.Value];

if ~app.NoLightButton.Value
    distance_matrix_I = NaN(frames_I, 7, n_I);
    orientationxy_matrix_I = NaN(frames_I, n_I);
    
    for i = 1:ninhibition
        dzlight_all_I{i} = [];
        dxyzlight_all_I{i} = [];
        orientationxylight_all_I{i} = [];
        
        dzlight_mean_I{i} = [];
        dxyzlight_mean_I{i} = [];
        orientationxylight_mean_I{i} = [];
        
        dzlight_var_I{i} = [];
        dxyzlight_var_I{i} = [];
        orientationxylight_var_I{i} = [];
        
        dz1_I{i} = [];
        dz2_I{i} = [];
        z1_I{i} = [];
        z2_I{i} = [];
    end
    for i = 1:n_I
        frames_temp = size(distance_all_I{i}, 1);
        distance_matrix_I(1:frames_temp, :, i) = distance_all_I{i};
        orientationxy_matrix_I(1:frames_temp, i) = orientationxy_all_I{i};
        
        for j = 1:ninhibition
            dzlight_all_I{j} = [dzlight_all_I{j}; distance_all_I{i}(laserstart_all(j, i):laserstop_all(j, i), 3)];
            dxyzlight_all_I{j} = [dxyzlight_all_I{j}; distance_all_I{i}(laserstart_all(j, i):laserstop_all(j, i), 7)];
            orientationxylight_all_I{j} = [orientationxylight_all_I{j}; orientationxy_all_I{i}(laserstart_all(j, i):laserstop_all(j, i))];
            
            dzlight_mean_I{j} = [dzlight_mean_I{j}; mean(distance_all_I{i}(laserstart_all(j, i):laserstop_all(j, i), 3), 'omitnan')];
            dxyzlight_mean_I{j} = [dxyzlight_mean_I{j}; mean(distance_all_I{i}(laserstart_all(j, i):laserstop_all(j, i), 7), 'omitnan')];
            orientationxylight_mean_I{j} = [orientationxylight_mean_I{j}; mean(orientationxy_all_I{i}(laserstart_all(j, i):laserstop_all(j, i)), 'omitnan')];
            
            dzlight_var_I{j} = [dzlight_var_I{j}; var(distance_all_I{i}(laserstart_all(j, i):laserstop_all(j, i), 3), 'omitnan')];
            dxyzlight_var_I{j} = [dxyzlight_var_I{j}; var(distance_all_I{i}(laserstart_all(j, i):laserstop_all(j, i), 7), 'omitnan')];
            orientationxylight_var_I{j} = [orientationxylight_var_I{j}; var(orientationxy_all_I{i}(laserstart_all(j, i):laserstop_all(j, i)), 'omitnan')];
            
            dz1_I{j} = [dz1_I{j}; dz1_all_I{i}(laserstart_all(j, i):laserstop_all(j, i))];
            dz2_I{j} = [dz2_I{j}; dz2_all_I{i}(laserstart_all(j, i):laserstop_all(j, i))];
            z1_I{j} = [z1_I{j}; z1_all_I{i}(laserstart_all(j, i):laserstop_all(j, i))];
            z2_I{j} = [z2_I{j}; z2_all_I{i}(laserstart_all(j, i):laserstop_all(j, i))];
        end
    end
    laserstart = mean(laserstart_all, 2)/FrameRate;
    laserstop = mean(laserstop_all, 2)/FrameRate;
    t_I = (1:frames_I)'/FrameRate;
    
    ylabel_text = {'\DeltaX (mm)', '\DeltaY (mm)', '\DeltaZ (mm)',...
        '\DeltaXY (mm)', '\DeltaXZ (mm)', '\DeltaYZ (mm)', '\DeltaXYZ (mm)'};
    title_text = {'\DeltaX', '\DeltaY', '\DeltaZ', '\DeltaXY', '\DeltaXZ', '\DeltaYZ', '\DeltaXYZ'};
    for i = [3 7]
        figure;
        plot_tj_multitrial(time_range, t_NI, squeeze(distance_matrix_NI(:, i, :)), t_I, squeeze(distance_matrix_I(:, i, :)),...
            laserstart, laserstop, ylabel_text{i}, title_text{i});
    end
    
    figure;
    plot_tj_multitrial(time_range, t_NI, orientationxy_matrix_NI, t_I, orientationxy_matrix_I,...
        laserstart, laserstop, ['Orientation (' char(176) ')'], 'REF XY Plane');
else
    ylabel_text = {'\DeltaX (mm)', '\DeltaY (mm)', '\DeltaZ (mm)',...
        '\DeltaXY (mm)', '\DeltaXZ (mm)', '\DeltaYZ (mm)', '\DeltaXYZ (mm)'};
    title_text = {'\DeltaX', '\DeltaY', '\DeltaZ', '\DeltaXY', '\DeltaXZ', '\DeltaYZ', '\DeltaXYZ'};
    for i = [3 7]
        figure;
        plot_tj_multitrial(time_range, t_NI, squeeze(distance_matrix_NI(:, i, :)), [], [],...
            [], [], ylabel_text{i}, title_text{i});
    end
    
    figure;
    plot_tj_multitrial(time_range, t_NI, orientationxy_matrix_NI, [], [],...
        laserstart, laserstop, ['Orientation (' char(176) ')'], 'REF XY Plane');
end

t = (-trange:1:trange)'/FrameRate;
time_range = [t(1) t(end)];
if ~app.NoLightButton.Value
    for i = 1:ninhibition
        % trajectories aligned to bite time
        figure;
        plot_tj_multitrial(time_range, t, dzbite_NI(:, 1:nbite_NI(i), i), t, dzbite_I(:, 1:nbite_I(i), i), 0, NaN, '\DeltaZ (mm)', ['Inhibition #' num2str(i)]);
        figure;
        plot_tj_multitrial(time_range, t, dxyzbite_NI(:, 1:nbite_NI(i), i), t, dxyzbite_I(:, 1:nbite_I(i), i), 0, NaN, '\DeltaXYZ (mm)', ['Inhibition #' num2str(i)]);
        figure;
        plot_tj_multitrial(time_range, t, orientationxybite_NI(:, 1:nbite_NI(i), i), t, orientationxybite_I(:, 1:nbite_I(i), i), 0, NaN, ['Orientation (' char(176) ')'], ['REF XY Plane (Inhibition #' num2str(i) ')']);
        
        % histograms
        if nbite_I(i) ~= 0
            histogram_multitrial(dzlight_all_NI{i}, dzbite_NI(trange+1, :, i), dzlight_all_I{i}, dzbite_I(trange+1, :, i), '\DeltaZ (mm)', ['\DeltaZ (Inhibition #' num2str(i) ')']);
            histogram_multitrial(dxyzlight_all_NI{i}, dxyzbite_NI(trange+1, :, i), dxyzlight_all_I{i}, dxyzbite_I(trange+1, :, i), '\DeltaXYZ (mm)', ['\DeltaXYZ (Inhibition #' num2str(i) ')']);
            histogram_multitrial(orientationxylight_all_NI{i}, orientationxybite_NI(trange+1, :, i), orientationxylight_all_I{i}, orientationxybite_I(trange+1, :, i),...
                ['Orientation (' char(176) ')'], ['REF XY Plane (Inhibition #' num2str(i) ')']);
        else
            histogram_multitrial(dzlight_all_NI{i}, dzbite_NI(trange+1, :, i), dzlight_all_I{i}, [], '\DeltaZ (mm)', ['\DeltaZ (Inhibition #' num2str(i) ')']);
            histogram_multitrial(dxyzlight_all_NI{i}, dxyzbite_NI(trange+1, :, i), dxyzlight_all_I{i}, [], '\DeltaXYZ (mm)', ['\DeltaXYZ (Inhibition #' num2str(i) ')']);
            histogram_multitrial(orientationxylight_all_NI{i}, orientationxybite_NI(trange+1, :, i), orientationxylight_all_I{i}, [],...
                ['Orientation (' char(176) ')'], ['REF XY Plane (Inhibition #' num2str(i) ')']);
        end
        
        histogram_multitrial(log(dzlight_var_NI{i}), [], log(dzlight_var_I{i}), [], 'Log \DeltaZ Variance (mm^{2})', ['\DeltaZ (Inhibition #' num2str(i) ')']);
        histogram_multitrial(log(dxyzlight_var_NI{i}), [], log(dxyzlight_var_I{i}), [], 'Log \DeltaXYZ Variance (mm^{2})', ['\DeltaXYZ (Inhibition #' num2str(i) ')']);
        histogram_multitrial(log(orientationxylight_var_NI{i}), [], log(orientationxylight_var_I{i}), [], ['Log Orientation Variance (' char(176) '^{2})'], ['REF XY Plane (Inhibition #' num2str(i) ')']);
        
        % correlations
%         z_threshold = zspeed_threshold/FrameRate;
%         dz1_NI_temp = dz1_NI{i}(~(abs(dz1_NI{i}) <= z_threshold & abs(dz2_NI{i}) <= z_threshold));
%         dz2_NI_temp = dz2_NI{i}(~(abs(dz1_NI{i}) <= z_threshold & abs(dz2_NI{i}) <= z_threshold));
%         dz1_I_temp = dz1_I{i}(~(abs(dz1_I{i}) <= z_threshold & abs(dz2_I{i}) <= z_threshold));
%         dz2_I_temp = dz2_I{i}(~(abs(dz1_I{i}) <= z_threshold & abs(dz2_I{i}) <= z_threshold));
        correlation_multitrial(dz1_NI{i}, dz2_NI{i}, dz1_I{i}, dz2_I{i}, [app.DropDown_2.Items{pointID1} ' \DeltaZ (mm)'], [app.DropDown_3.Items{pointID2} ' \DeltaZ (mm)'], ['\DeltaZ Correlation (Inhibition #' num2str(i) ')'], 'equal');
        correlation_multitrial(z1_NI{i}, z2_NI{i}, z1_I{i}, z2_I{i}, [app.DropDown_2.Items{pointID1} ' Z (mm)'], [app.DropDown_3.Items{pointID2} ' Z (mm)'], ['Z Correlation (Inhibition #' num2str(i) ')'], 'equal');
        
        result.dzlight_mean_NI(i) = mean(dzlight_mean_NI{i}, 'omitnan');
        result.dxyzlight_mean_NI(i) = mean(dxyzlight_mean_NI{i}, 'omitnan');
        result.orientationxylight_mean_NI(i) = mean(orientationxylight_mean_NI{i}, 'omitnan');
        
        result.dzlight_var_NI(i) = mean(dzlight_var_NI{i}, 'omitnan');
        result.dxyzlight_var_NI(i) = mean(dxyzlight_var_NI{i}, 'omitnan');
        result.orientationxylight_var_NI(i) = mean(orientationxylight_var_NI{i}, 'omitnan');
        
        result.dzlight_mean_I(i) = mean(dzlight_mean_I{i}, 'omitnan');
        result.dxyzlight_mean_I(i) = mean(dxyzlight_mean_I{i}, 'omitnan');
        result.orientationxylight_mean_I(i) = mean(orientationxylight_mean_I{i}, 'omitnan');
        
        result.dzlight_var_I(i) = mean(dzlight_var_I{i}, 'omitnan');
        result.dxyzlight_var_I(i) = mean(dxyzlight_var_I{i}, 'omitnan');
        result.orientationxylight_var_I(i) = mean(orientationxylight_var_I{i}, 'omitnan');
        
        result.dzbite_NI{i} = dzbite_NI(:, 1:nbite_NI(i), i);
        result.dxyzbite_NI{i} = dxyzbite_NI(:, 1:nbite_NI(i), i);
        result.orientationxybite_NI{i} = orientationxybite_NI(:, 1:nbite_NI(i), i);
        
        result.dzbite_I{i} = dzbite_I(:, 1:nbite_I(i), i);
        result.dxyzbite_I{i} = dxyzbite_I(:, 1:nbite_I(i), i);
        result.orientationxybite_I{i} = orientationxybite_I(:, 1:nbite_I(i), i);
    end
else
    for i = 1:ninhibition
        % trajectories aligned to bite time
        figure;
        plot_tj_multitrial(time_range, t, dzbite_NI(:, 1:nbite_NI(i), i), [], [], 0, NaN, '\DeltaZ (mm)', ['Inhibition #' num2str(i)]);
        figure;
        plot_tj_multitrial(time_range, t, dxyzbite_NI(:, 1:nbite_NI(i), i), [], [], 0, NaN, '\DeltaXYZ (mm)', ['Inhibition #' num2str(i)]);
        figure;
        plot_tj_multitrial(time_range, t, orientationxybite_NI(:, 1:nbite_NI(i), i), [], [], 0, NaN, ['Orientation (' char(176) ')'], ['REF XY Plane (Inhibition #' num2str(i) ')']);
        
        % histograms
        histogram_multitrial(dzlight_all_NI{i}, dzbite_NI(trange+1, :, i),[], [], '\DeltaZ (mm)', ['\DeltaZ (Inhibition #' num2str(i) ')']);
        histogram_multitrial(dxyzlight_all_NI{i}, dxyzbite_NI(trange+1, :, i), [], [], '\DeltaXYZ (mm)', ['\DeltaXYZ (Inhibition #' num2str(i) ')']);
        histogram_multitrial(orientationxylight_all_NI{i}, orientationxybite_NI(trange+1, :, i), [], [], ['Orientation (' char(176) ')'], ['REF XY Plane (Inhibition #' num2str(i) ')']);
        
        histogram_multitrial(log(dzlight_var_NI{i}), [], [], [], 'Log \DeltaZ Variance (mm^{2})', ['\DeltaZ (Inhibition #' num2str(i) ')']);
        histogram_multitrial(log(dxyzlight_var_NI{i}), [], [], [], 'Log \DeltaXYZ Variance (mm^{2})', ['\DeltaXYZ (Inhibition #' num2str(i) ')']);
        histogram_multitrial(log(orientationxylight_var_NI{i}), [], [], [], ['Log Orientation Variance (' char(176) '^{2})'], ['REF XY Plane (Inhibition #' num2str(i) ')']);
        
        correlation_multitrial(dz1_NI{i}, dz2_NI{i}, [], [], [app.DropDown_2.Items{pointID1} ' \DeltaZ (mm)'], [app.DropDown_3.Items{pointID2} ' \DeltaZ (mm)'], ['\DeltaZ Correlation (Inhibition #' num2str(i) ')'], 'equal');
        correlation_multitrial(z1_NI{i}, z2_NI{i}, [], [], [app.DropDown_2.Items{pointID1} ' Z (mm)'], [app.DropDown_3.Items{pointID2} ' Z (mm)'], ['Z Correlation (Inhibition #' num2str(i) ')'], 'equal');
        
        result.dzlight_mean_NI(i) = mean(dzlight_mean_NI{i}, 'omitnan');
        result.dxyzlight_mean_NI(i) = mean(dxyzlight_mean_NI{i}, 'omitnan');
        result.orientationxylight_mean_NI(i) = mean(orientationxylight_mean_NI{i}, 'omitnan');
        
        result.dzlight_var_NI(i) = mean(dzlight_var_NI{i}, 'omitnan');
        result.dxyzlight_var_NI(i) = mean(dxyzlight_var_NI{i}, 'omitnan');
        result.orientationxylight_var_NI(i) = mean(orientationxylight_var_NI{i}, 'omitnan');
        
        result.dzbite_NI{i} = dzbite_NI(:, 1:nbite_NI(i), i);
        result.dxyzbite_NI{i} = dxyzbite_NI(:, 1:nbite_NI(i), i);
        result.orientationxybite_NI{i} = orientationxybite_NI(:, 1:nbite_NI(i), i);
    end
end

[file, path] = uiputfile('AverageDistance.mat', 'Save result');
if file ~= 0
    save([path file], 'result');
end