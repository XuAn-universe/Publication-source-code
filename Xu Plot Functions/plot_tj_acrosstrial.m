function plot_tj_acrosstrial(t, yall, lcolor, pcolor, ylabel_text)
hold on;
meany = mean(yall, 2, 'omitnan');
semy = std(yall, 0, 2, 'omitnan')/sqrt(size(yall, 2));
patch([t; t(end:-1:1)], [meany+semy; meany(end:-1:1)-semy(end:-1:1)], pcolor, 'EdgeColor', 'none', 'FaceAlpha', 0.2);
plot(t, meany, '-', 'Color', lcolor);
xlabel('Time (s)');
ylabel(ylabel_text);
box off;
set(gca, 'TickDir', 'out', 'FontSize', 12, 'ButtonDownFcn', @extract_figure);