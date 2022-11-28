function singlebar_plot(data, x, XTickLabel_text, YLabel_text)
hold on;
bar(x, mean(data), 1, 'FaceColor', [0.5 0.5 0.5]);
plot(x, data, 'o', 'MarkerFaceColor', [1 1 1], 'MarkerEdgeColor', [0 0 0]);
errorbar(x, mean(data), std(data)/sqrt(numel(data)), 'k', 'LineStyle', 'none', 'CapSize', 15);
set(gca, 'XTick', x, 'XTickLabel', XTickLabel_text, 'XLim', [0 x+1], 'TickLength', [0 0], 'FontSize', 12);
ylabel(YLabel_text);
box off;