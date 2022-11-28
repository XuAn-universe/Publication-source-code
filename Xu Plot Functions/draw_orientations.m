function draw_orientations(orientation_NI, orientation_I, edges, title_text)
if ~isempty(orientation_NI)
    histogram(orientation_NI, edges, 'FaceColor', [0 0 0], 'FaceAlpha', 1, 'Normalization', 'probability');
end
hold on;
if ~isempty(orientation_I)
    histogram(orientation_I, edges, 'FaceColor', 'none', 'EdgeColor', [0 1 0], 'FaceAlpha', 1, 'Normalization', 'probability');
end
xlim([edges(1) edges(end)]);
set(gca, 'TickDir', 'out', 'FontSize', 12, 'ButtonDownFcn', @extract_figure);
box off;
xlabel(['Orientation (' char(176) ')']);
ylabel('Probability');
title(title_text);
end