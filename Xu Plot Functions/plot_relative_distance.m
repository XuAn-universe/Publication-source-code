function plot_relative_distance(app, Exp_Path, FrameRate, pre, post)
% FrameRate = 120; %
value = app.VideosListBox.Value;

if app.TrackingDataCheckBox.Value
    label_table = table2array(app.UITable.Data);
else
    helpdlg('''Tracking Data'' needs to be checked');
    return;
end

onframes = zeros(1, numel(value));
pawl2nose_all = cell(1, numel(value));
pawr2nose_all = cell(1, numel(value));
pawl2pawr_all = cell(1, numel(value));

for i = 1:numel(value)
    csvdata = load([Exp_Path '\' app.VideosListBox.Items{value(i)}(1:3) 'gpio' app.VideosListBox.Items{value(i)}(10:end-4) '.csv']);
    startID = find(csvdata(:, 1) == 1, 1, 'first');
    stopID = find(csvdata(:, 1) == 1, 1, 'last');
    onframes(i) = stopID-startID+1;
    
    camtrack = find_tracking_filename(Exp_Path, app.VideosListBox.Items{value(i)}(1:end-4));
    if ~isempty(camtrack)
        [num_track, ~, ~] = xlsread(camtrack);
    end
    noseID = 3;
    nose_tj = num_track(:, noseID*3-1:noseID*3);
    nose_tj(num_track(:, noseID*3+1) < label_table(noseID, 1)) = nan;
    nose_tj_ROI = nose_tj(startID-pre:stopID+post, :);
    pawlID = 4;
    pawl_tj = num_track(:, pawlID*3-1:pawlID*3);
    pawl_tj(num_track(:, pawlID*3+1) < label_table(pawlID, 1)) = nan;
    pawl_tj_ROI = pawl_tj(startID-pre:stopID+post, :);
    pawrID = 5;
    pawr_tj = num_track(:, pawrID*3-1:pawrID*3);
    pawr_tj(num_track(:, pawrID*3+1) < label_table(pawrID, 1)) = nan;
    pawr_tj_ROI = pawr_tj(startID-pre:stopID+post, :);
    
    pawl2nose = vecnorm(pawl_tj_ROI-nose_tj_ROI, 2, 2);
    pawr2nose = vecnorm(pawr_tj_ROI-nose_tj_ROI, 2, 2);
    pawl2pawr = vecnorm(pawl_tj_ROI-pawr_tj_ROI, 2, 2);
    pawl2nose_all{i} = pawl2nose;
    pawr2nose_all{i} = pawr2nose;
    pawl2pawr_all{i} = pawl2pawr;
end

pawl2nose = nan(pre+max(onframes)+post, numel(value));
pawr2nose = nan(pre+max(onframes)+post, numel(value));
pawl2pawr = nan(pre+max(onframes)+post, numel(value));
t = ((1:pre+max(onframes)+post)-pre-1)'/FrameRate;
figure;
hold on;
for i = 1:numel(value)
    plot(t(1:numel(pawl2nose_all{i})), pawl2nose_all{i}, '-', 'Color', [0 187 255]/255, 'LineWidth', 1);
    plot(t(1:numel(pawr2nose_all{i})), pawr2nose_all{i}, '-', 'Color', [0 255 0]/255, 'LineWidth', 1);
    
    pawl2nose(1:numel(pawl2nose_all{i}), i) = pawl2nose_all{i};
    pawr2nose(1:numel(pawr2nose_all{i}), i) = pawr2nose_all{i};
end
hp(1) = plot(t, mean(pawl2nose, 2, 'omitnan'), '-', 'Color', [0 0 255]/255, 'LineWidth', 2);
hp(2) = plot(t, mean(pawr2nose, 2, 'omitnan'), '-', 'Color', [0 130 0]/255, 'LineWidth', 2);
yl = ylim;
patch([0 0 mean(onframes)/FrameRate mean(onframes)/FrameRate], [0 yl(2) yl(2) 0], [0.9 0.9 0.9], 'EdgeColor', 'none');
set(gca, 'TickDir', 'out', 'FontSize', 12);
hc = get(gca, 'Children');
set(gca, 'Children', [hc(2:end); hc(1)]);
xlabel('Time (s)');
ylabel('Hand-to-nose distance (pixel)');
legend(hp, {'Contra', 'Ipsi'}, 'FontSize', 12);
legend('boxoff');

figure;
hold on;
for i = 1:numel(value)
    plot(t(1:numel(pawl2pawr_all{i})), pawl2pawr_all{i}, '-', 'Color', [0.8 0.8 0.8], 'LineWidth', 1);
    
    pawl2pawr(1:numel(pawl2pawr_all{i}), i) = pawl2pawr_all{i};
end
plot(t, mean(pawl2pawr, 2, 'omitnan'), '-', 'Color', [0 0 0]/255, 'LineWidth', 2);
yl = ylim;
patch([0 0 mean(onframes)/FrameRate mean(onframes)/FrameRate], [0 yl(2) yl(2) 0], [0.9 0.9 0.9], 'EdgeColor', 'none');
set(gca, 'TickDir', 'out', 'FontSize', 12);
hc = get(gca, 'Children');
set(gca, 'Children', [hc(2:end); hc(1)]);
xlabel('Time (s)');
ylabel('Hand-to-hand distance (pixel)');

figure;
hold on;
hp(1) = plot_tj_MeanSEM(t, pawl2nose, [0 0 255]/255, [0 0 255]/255, '', '', '');
hp(2) = plot_tj_MeanSEM(t, pawr2nose, [0 130 0]/255, [0 130 0]/255, 'Time (s)', 'Hand-to-nose distance (pixel)', '');
yl = ylim;
patch([0 0 mean(onframes)/FrameRate mean(onframes)/FrameRate], [0 yl(2) yl(2) 0], [0.9 0.9 0.9], 'EdgeColor', 'none');
hc = get(gca, 'Children');
set(gca, 'Children', [hc(2:end); hc(1)]);
legend(hp, {'Contra', 'Ipsi'}, 'FontSize', 12);
legend('boxoff');

figure;
hold on;
plot_tj_MeanSEM(t, pawl2pawr, [255 153 0]/255, [255 153 0]/255, 'Time (s)', 'Hand-to-hand distance (pixel)', '');
yl = ylim;
patch([0 0 mean(onframes)/FrameRate mean(onframes)/FrameRate], [0 yl(2) yl(2) 0], [0.9 0.9 0.9], 'EdgeColor', 'none');
hc = get(gca, 'Children');
set(gca, 'Children', [hc(2:end); hc(1)]);

% set(gca, 'XLim', [-0.3 0.8], 'XTick', [-0.3 0 0.5 0.8]);
