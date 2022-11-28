function hp = plot_tj_MeanSEM(x, y, patch_color, mean_color, xlabel_text, ylabel_text, title_text)
% for y, column is repetition
if size(x, 2) > 1
    x = x';
end
ntrial = size(y, 2);
y_mean = mean(y, 2, 'omitnan');
y_sem = std(y, 0, 2, 'omitnan')/sqrt(ntrial);
hold on;
patch([x; x(end:-1:1)], [y_mean+y_sem; y_mean(end:-1:1)-y_sem(end:-1:1)], patch_color, 'EdgeColor', 'none', 'FaceAlpha', 0.2);
hp = plot(x, y_mean, '-', 'Color', mean_color, 'LineWidth', 1);
[xl(1), xl(2)] = bounds(x);
xlim(xl);
% [yl(1), yl(2)] = bounds([y_mean+y_sem; y_mean-y_sem]);
% ylim(yl);
if ~isempty(xlabel_text)
    xlabel(xlabel_text);
end
if ~isempty(ylabel_text)
    ylabel(ylabel_text);
end
if ~isempty(title_text)
    title(title_text);
end
box off;
set(gca, 'TickDir', 'out', 'FontSize', 12, 'ButtonDownFcn', @extract_figure);