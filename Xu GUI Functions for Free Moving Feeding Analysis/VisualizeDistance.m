function VisualizeDistance(app, Exp_Path, FrameRate, nmedian, pawsegment)
trial = app.TrialsListBox.Value;
if numel(trial) ~= 1
    helpdlg('You can only choose one trial');
    return;
end

trange = app.bitesEditField.Value; % time range for bite alignment
trange = round(trange*FrameRate);

nbite = 0;
zbite = [];
orientationxybite = [];

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
% process trajectories
if app.TrackingDataCheckBox.Value
    label_table = table2array(app.UITable.Data);
else
    helpdlg('''Tracking Data'' needs to be checked');
    return;
end
pointID1 = app.DropDown_2.Value;
pointID2 = app.DropDown_3.Value;

[~, ~, ~, ~, ~, x1, y1, z1, ~, ~, laserstart, laserstop] = trajectory_postprocessing(pointID1, Exp_Path, trial,...
    label_table, nmedian, FrameRate);
[~, ~, ~, ~, ~, x2, y2, z2, ~, ~, ~, ~] = trajectory_postprocessing(pointID2, Exp_Path, trial, label_table, nmedian, FrameRate);
dx = x1(:, 1)-x2(:, 1);
dx_filtered = medfilt1(dx, nmedian, 'omitnan', 'truncate');
dy = y1(:, 1)-y2(:, 1);
dy_filtered = medfilt1(dy, nmedian, 'omitnan', 'truncate');
dz = z1(:, 1)-z2(:, 1);
dz_filtered = medfilt1(dz, nmedian, 'omitnan', 'truncate');
dxy = sqrt(dx.^2+dy.^2);
dxy_filtered = medfilt1(dxy, nmedian, 'omitnan', 'truncate');
dxz = sqrt(dx.^2+dz.^2);
dxz_filtered = medfilt1(dxz, nmedian, 'omitnan', 'truncate');
dyz = sqrt(dy.^2+dz.^2);
dyz_filtered = medfilt1(dyz, nmedian, 'omitnan', 'truncate');
dxyz = sqrt(dx.^2+dy.^2+dz.^2);
dxyz_filtered = medfilt1(dxyz, nmedian, 'omitnan', 'truncate');
distances = [dx dy dz dxy dxz dyz dxyz];
distances_filtered = [dx_filtered dy_filtered dz_filtered dxy_filtered dxz_filtered dyz_filtered dxyz_filtered];
orientations = [atan2d(dy, dx) atan2d(dz, dx) atan2d(dz, dy)];
orientations_filtered = medfilt1(orientations, nmedian, 'omitnan', 'truncate');
orientationxy = atand(dz./(sqrt(dx.^2+dy.^2)));
orientationxy_filtered = medfilt1(orientationxy, nmedian, 'omitnan', 'truncate');
t = (1:size(x1, 1))'/FrameRate;
time_range = [app.timerangesEditField.Value app.toEditField.Value];
if app.NoLightButton.Value || ~app.OptSessionCheckBox.Value
    laserstart = [];
    laserstop = [];
else
    laserstart = laserstart/FrameRate;
    laserstop = laserstop/FrameRate;
end

if ~isempty(bite_timestamps)
%     figure;
%     hold on;
    for i = 1:numel(bite_timestamps)
        timestamp = bite_timestamps(i);
        [~, id] = min(abs(t-timestamp));
        if id-trange >= 1
            nbite = nbite+1;
            zbite(:, nbite) = z2(id-trange:id+trange, 1);
            orientationxybite(:, nbite) = orientationxy(id-trange:id+trange);
            
            if 0
                figure;
                colorfulplot(orientationxybite(:, nbite)-orientationxybite(trange+1, nbite), zbite(:, nbite)-zbite(trange+1, nbite), ['Orientation (' char(176) ')'], 'Z (mm)', [], 'jet');
                %             colorfulplot(orientationxybite(:, nbite)-min(orientationxybite(:, nbite)), zbite(:, nbite)-max(zbite(:, nbite)));
                hold on;
                plot(0, 0, '.k', 'MarkerSize', 20);
            end
        end
    end
%     plot(0, 0, '.k', 'MarkerSize', 20);
end

% raw and filtered X, Y, Z, XY, XZ, YZ, XYZ distances
label_distance = {'\DeltaX', '\DeltaY', '\DeltaZ', '\DeltaXY', '\DeltaXZ', '\DeltaYZ', '\DeltaXYZ'};
label_distance_filtered = {'\DeltaX filtered', '\DeltaY filtered', '\DeltaZ filtered', '\DeltaXY filtered',...
    '\DeltaXZ filtered', '\DeltaYZ filtered', '\DeltaXYZ filtered'};
figure;
for i = 1:7
    subplot(3, 3, i);
    plot_tj_singletrial(t, [distances(:, i) distances_filtered(:, i)], time_range, laserstart, laserstop,...
        [0 1 0], {'r', 'k'}, 'distance (mm)', {label_distance{i}, label_distance_filtered{i}});
end

% raw X, Y, Z, XY, XZ, YZ, XYZ distances with bite events
label_distance = {'\DeltaX (mm)', '\DeltaY (mm)', '\DeltaZ (mm)',...
    '\DeltaXY (mm)', '\DeltaXZ (mm)', '\DeltaYZ (mm)', '\DeltaXYZ (mm)'};
figure;
for i = 1:7
    subplot(3, 3, i);
    plot_tj_singletrial(t, distances(:, i), time_range, laserstart, laserstop, [0 1 0], {'k'}, label_distance{i}, '');
    if ~isempty(bite_timestamps)
        plot_biteevents(bite_timestamps, bite_amplitudes);
    end
end

% raw and filtered orientation with and without bite events
label_orientation = {'XY Plane', 'XZ Plane', 'YZ Plane'};
label_orientation_filtered = {'XY Plane filtered', 'XZ Plane filtered', 'YZ Plane filtered'};
figure;
for i = 1:3
    subplot(2, 3, i);
    plot_tj_singletrial(t, [orientations(:, i) orientations_filtered(:, i)], time_range, laserstart, laserstop,...
        [0 1 0], {'r', 'k'}, ['Orientation (' char(176) ')'], {label_orientation{i}, label_orientation_filtered{i}});
end
for i = 1:3
    subplot(2, 3, i+3);
    plot_tj_singletrial(t, orientations(:, i), time_range, laserstart, laserstop,...
        [0 1 0], {'k'}, ['Orientation (' char(176) ')'], label_orientation{i});
    if ~isempty(bite_timestamps)
        plot_biteevents(bite_timestamps, bite_amplitudes);
    end
    title(label_orientation{i});
end

figure;
subplot(1, 2, 1);
plot_tj_singletrial(t, [orientationxy orientationxy_filtered], time_range, laserstart, laserstop,...
    [0 1 0], {'r', 'k'}, ['Orientation (' char(176) ')'], {'XY', 'XY filtered'});
subplot(1, 2, 2);
plot_tj_singletrial(t, orientationxy, time_range, laserstart, laserstop,...
    [0 1 0], {'k'}, ['Orientation (' char(176) ')'], 'XY');
if ~isempty(bite_timestamps)
    plot_biteevents(bite_timestamps, bite_amplitudes);
end

% polar plot for relative locations
title_text = {'XY Plane', 'XZ Plane', 'YZ Plane'};
figure;
for i = 1:3
    subplot(2, 3, i);
    polarplot_tj_singletrial(orientations(:, i), distances(:, i+3), t, time_range, laserstart, laserstop, bite_timestamps);
    title(title_text{i});
end
title_text = {'XY Plane filtered', 'XZ Plane filtered', 'YZ Plane filtered'};
for i = 1:3
    subplot(2, 3, i+3);
    polarplot_tj_singletrial(orientations_filtered(:, i), distances_filtered(:, i+3),...
        t, time_range, laserstart, laserstop, bite_timestamps);
    title(title_text{i});
end

return; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% trajectory in 3D
x1 = x1(t >= time_range(1) & t <= time_range(2), :);
y1 = y1(t >= time_range(1) & t <= time_range(2), :);
z1 = z1(t >= time_range(1) & t <= time_range(2), :);
x2 = x2(t >= time_range(1) & t <= time_range(2), :);
y2 = y2(t >= time_range(1) & t <= time_range(2), :);
z2 = z2(t >= time_range(1) & t <= time_range(2), :);
t = t(t >= time_range(1) & t <= time_range(2));
plot_lines(x1, y1, z1, x2, y2, z2, t, laserstart, laserstop);

% trajectory in 3D
x1 = downsample(x1, pawsegment);
y1 = downsample(y1, pawsegment);
z1 = downsample(z1, pawsegment);
x2 = downsample(x2, pawsegment);
y2 = downsample(y2, pawsegment);
z2 = downsample(z2, pawsegment);
t = downsample(t, pawsegment);
plot_lines(x1, y1, z1, x2, y2, z2, t, laserstart, laserstop);

    function plot_lines(x1, y1, z1, x2, y2, z2, t, laserstart, laserstop)
        figure;
        haxis = zeros(1, 4);
        xlabel_text = {'X (mm)', 'X (mm)', 'Y (mm)', 'X (mm)'};
        ylabel_text = {'Y (mm)', 'Z (mm)', 'Z (mm)', 'Y (mm)'};
        for ii = 1:4
            haxis(ii) = subplot(2, 2, ii);
            hold on;
            xlabel(xlabel_text{ii});
            ylabel(ylabel_text{ii});
            if ii == 4
                zlabel('Z (mm)');
            end
            box off;
        end
        for ii = 1:numel(t)
            x1_temp = x1(ii, 1);
            y1_temp = y1(ii, 1);
            z1_temp = z1(ii, 1);
            x2_temp = x2(ii, 1);
            y2_temp = y2(ii, 1);
            z2_temp = z2(ii, 1);
            if ~isempty(laserstart)
                for j = 1:numel(laserstart)
                    if t(ii) >= laserstart(j) && t(ii) <= laserstop(j)
                        lcolor = [0 1 0];
                    else
                        lcolor = [0 0 0];
                    end
                end
            else
                lcolor = [0 0 0];
            end
            line([x1_temp; x2_temp], [y1_temp; y2_temp], 'Parent', haxis(1), 'Color', lcolor, 'LineWidth', 0.5);
            plot(x1_temp, y1_temp, 'Parent', haxis(1), 'Marker', 'o', 'MarkerEdgeColor', lcolor, 'MarkerFaceColor', 'none', 'LineStyle', 'none');
            plot(x2_temp, y2_temp, 'Parent', haxis(1), 'Marker', 's', 'MarkerEdgeColor', lcolor, 'MarkerFaceColor', 'none', 'LineStyle', 'none');
            line([x1_temp; x2_temp], [z1_temp; z2_temp], 'Parent', haxis(2), 'Color', lcolor, 'LineWidth', 0.5);
            plot(x1_temp, z1_temp, 'Parent', haxis(2), 'Marker', 'o', 'MarkerEdgeColor', lcolor, 'MarkerFaceColor', 'none', 'LineStyle', 'none');
            plot(x2_temp, z2_temp, 'Parent', haxis(2), 'Marker', 's', 'MarkerEdgeColor', lcolor, 'MarkerFaceColor', 'none', 'LineStyle', 'none');
            line([y1_temp; y2_temp], [z1_temp; z2_temp], 'Parent', haxis(3), 'Color', lcolor, 'LineWidth', 0.5);
            plot(y1_temp, z1_temp, 'Parent', haxis(3), 'Marker', 'o', 'MarkerEdgeColor', lcolor, 'MarkerFaceColor', 'none', 'LineStyle', 'none');
            plot(y2_temp, z2_temp, 'Parent', haxis(3), 'Marker', 's', 'MarkerEdgeColor', lcolor, 'MarkerFaceColor', 'none', 'LineStyle', 'none');
            line([x1_temp; x2_temp], [y1_temp; y2_temp], [z1_temp; z2_temp], 'Parent', haxis(4), 'Color', lcolor, 'LineWidth', 0.5);
            plot3(x1_temp, y1_temp, z1_temp, 'Parent', haxis(4), 'Marker', 'o', 'MarkerEdgeColor', lcolor, 'MarkerFaceColor', 'none', 'LineStyle', 'none');
            plot3(x2_temp, y2_temp, z2_temp, 'Parent', haxis(4), 'Marker', 's', 'MarkerEdgeColor', lcolor, 'MarkerFaceColor', 'none', 'LineStyle', 'none');
        end
        set(gca, 'TickLength', [0 0], 'FontSize', 12);
    end
end