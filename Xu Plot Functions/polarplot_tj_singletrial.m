function polarplot_tj_singletrial(orientations, distances, t, time_range, laserstart, laserstop, bite_timestamps)
orientations = orientations(t >= time_range(1) & t <= time_range(2));
distances = distances(t >= time_range(1) & t <= time_range(2));
t = t(t >= time_range(1) & t <= time_range(2));

polarplot(orientations/180*pi, distances, 'LineStyle', 'none', 'Marker', 'o', 'MarkerEdgeColor', [0 0 0]);
hold on;
if ~isempty(laserstart)
    orientations_I = [];
    distances_I = [];
    for i = 1:numel(laserstart)
        orientations_I = [orientations_I; orientations(t >= laserstart(i) & t <= laserstop(i))];
        distances_I = [distances_I; distances(t >= laserstart(i) & t <= laserstop(i))];
    end
    polarplot(orientations_I/180*pi, distances_I, 'LineStyle', 'none', 'Marker', 'o', 'MarkerEdgeColor', [0 1 0]);
end
if ~isempty(bite_timestamps)
    id_all = [];
    for i = 1:numel(bite_timestamps)
        timestamp = bite_timestamps(i);
        if timestamp >= t(1) && timestamp <= t(end)
            [~, id] = min(abs(t-timestamp));
            id_all = [id_all; id];
        end
    end
    polarplot(orientations(id_all)/180*pi, distances(id_all), 'LineStyle', 'none',...
        'Marker', 'o', 'MarkerFaceColor', [1 0 0], 'MarkerEdgeColor', 'none');
end
set(gca, 'TickLength', [0 0], 'FontSize', 12);
end