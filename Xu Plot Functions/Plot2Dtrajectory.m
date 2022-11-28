function Plot2Dtrajectory(csvdata, num_track, TableData, pre, post)
startID = find(csvdata(:, 1) == 1, 1, 'first');
stopID = find(csvdata(:, 1) == 1, 1, 'last');
nframes = pre+post+stopID-startID+1;

gap = 1;
linewidth = 2.5;
figure;
axis image;
hold on;
for i = 1:nframes
    if mod(i-1, gap) == 0
        line([num_track(startID-pre+i-1, 1*3-1) num_track(startID-pre+i-1, 2*3-1)],...
            [num_track(startID-pre+i-1, 1*3) num_track(startID-pre+i-1, 2*3)], 'Color', repmat(1-i/nframes, 1, 3), 'LineWidth', linewidth);
        line([num_track(startID-pre+i-1, 2*3-1) num_track(startID-pre+i-1, 3*3-1)],...
            [num_track(startID-pre+i-1, 2*3) num_track(startID-pre+i-1, 3*3)], 'Color', repmat(1-i/nframes, 1, 3), 'LineWidth', linewidth);
        line([num_track(startID-pre+i-1, 1*3-1) num_track(startID-pre+i-1, 3*3-1)],...
            [num_track(startID-pre+i-1, 1*3) num_track(startID-pre+i-1, 3*3)], 'Color', repmat(1-i/nframes, 1, 3), 'LineWidth', linewidth);
    end
    
    if i ~= nframes
        line([num_track(startID-pre+i-1, 4*3-1) num_track(startID-pre+i, 4*3-1)],...
            [num_track(startID-pre+i-1, 4*3) num_track(startID-pre+i, 4*3)], 'Color', hsv2rgb([0.68 i/nframes 1]), 'LineWidth', linewidth);
        line([num_track(startID-pre+i-1, 5*3-1) num_track(startID-pre+i, 5*3-1)],...
            [num_track(startID-pre+i-1, 5*3) num_track(startID-pre+i, 5*3)], 'Color', hsv2rgb([0.37 i/nframes 1]), 'LineWidth', linewidth);
    end
end
plot(num_track(stopID+post, 1*3-1), num_track(stopID+post, 1*3), 's', 'MarkerSize', 11, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k');
plot(num_track(stopID+post, 2*3-1), num_track(stopID+post, 2*3), 's', 'MarkerSize', 11, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k');
plot(num_track(stopID+post, 3*3-1), num_track(stopID+post, 3*3), 's', 'MarkerSize', 11, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k');
plot(num_track(startID-pre, 1*3-1), num_track(startID-pre, 1*3), 'o', 'MarkerSize', 8, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k');
plot(num_track(startID-pre, 2*3-1), num_track(startID-pre, 2*3), 'o', 'MarkerSize', 8, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k');
plot(num_track(startID-pre, 3*3-1), num_track(startID-pre, 3*3), 'o', 'MarkerSize', 8, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k');
plot(num_track(stopID+post, 4*3-1), num_track(stopID+post, 4*3), 's', 'MarkerSize', 11, 'MarkerFaceColor', hsv2rgb([0.68 1 1]), 'MarkerEdgeColor', hsv2rgb([0.68 1 1]));
plot(num_track(stopID+post, 5*3-1), num_track(stopID+post, 5*3), 's', 'MarkerSize', 11, 'MarkerFaceColor', hsv2rgb([0.37 1 1]), 'MarkerEdgeColor', hsv2rgb([0.37 1 1]));
plot(num_track(startID-pre, 4*3-1), num_track(startID-pre, 4*3), 'o', 'MarkerSize', 8, 'MarkerFaceColor', hsv2rgb([0.68 1 1]), 'MarkerEdgeColor', hsv2rgb([0.68 1 1]));
plot(num_track(startID-pre, 5*3-1), num_track(startID-pre, 5*3), 'o', 'MarkerSize', 8, 'MarkerFaceColor', hsv2rgb([0.37 1 1]), 'MarkerEdgeColor', hsv2rgb([0.37 1 1]));
set(gca, 'YDir', 'reverse');