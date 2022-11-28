function plothistogram(data, edges, ecolor, fcolor, xlabel_text)
histogram(data, edges, 'FaceColor', fcolor, 'EdgeColor', ecolor, 'Normalization', 'probability');
set(gca, 'TickDir', 'out', 'FontSize', 12, 'ButtonDownFcn', @extract_figure);
box off;
xlabel(xlabel_text);
ylabel('Probability');