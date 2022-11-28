function plot_tj_multitrial(time_range, t_NI, Y_NI, t_I, Y_I, laserstart, laserstop, ylabel_text, title_text)
hold on;
Y_NI = Y_NI(t_NI >= time_range(1) & t_NI <= time_range(2), :);
t_NI = t_NI(t_NI >= time_range(1) & t_NI <= time_range(2));
mean_Y_NI = mean(Y_NI, 2, 'omitnan');
sem_Y_NI = std(Y_NI, 0, 2, 'omitnan')/sqrt(size(Y_NI, 2));
Y_range = [mean_Y_NI+sem_Y_NI; mean_Y_NI-sem_Y_NI];
patch([t_NI; t_NI(end:-1:1)], [mean_Y_NI+sem_Y_NI; mean_Y_NI(end:-1:1)-sem_Y_NI(end:-1:1)], [0 0 0], 'EdgeColor', 'none', 'FaceAlpha', 0.2);

if ~isempty(Y_I)
    Y_I = Y_I(t_I >= time_range(1) & t_I <= time_range(2), :);
    t_I = t_I(t_I >= time_range(1) & t_I <= time_range(2));
    mean_Y_I = mean(Y_I, 2, 'omitnan');
    sem_Y_I = std(Y_I, 0, 2, 'omitnan')/sqrt(size(Y_I, 2));
    Y_range = [mean_Y_NI+sem_Y_NI; mean_Y_NI-sem_Y_NI; mean_Y_I+sem_Y_I; mean_Y_I-sem_Y_I];
    patch([t_I; t_I(end:-1:1)], [mean_Y_I+sem_Y_I; mean_Y_I(end:-1:1)-sem_Y_I(end:-1:1)], [0 1 0], 'EdgeColor', 'none', 'FaceAlpha', 0.2);
    plot(t_I, mean_Y_I, '-', 'Color', [0 1 0]);
end

plot(t_NI, mean_Y_NI, '-', 'Color', [0 0 0]);

if ~isempty(laserstart)
    for ii = 1:numel(laserstart)
        line([laserstart(ii); laserstart(ii)], [min(Y_range) max(Y_range)], 'Color', [0 0 0], 'LineStyle', '--', 'LineWidth', 1);
        line([laserstop(ii); laserstop(ii)], [min(Y_range) max(Y_range)], 'Color', [0 0 0], 'LineStyle', '--', 'LineWidth', 1);
    end
end
xlim(time_range);
ylim([min(Y_range) max(Y_range)]);
xlabel('Time (s)');
ylabel(ylabel_text);
title(title_text);
box off;
set(gca, 'TickDir', 'out', 'FontSize', 12, 'ButtonDownFcn', @extract_figure);
end