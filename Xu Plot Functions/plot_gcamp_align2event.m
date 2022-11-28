function plot_gcamp_align2event(t, zscore_random, zscore, title_text)
if ~isempty(zscore_random) && ~isempty(zscore)
    nevent = size(zscore, 2);
    zscore_random_mean = mean(zscore_random, 2, 'omitnan');
    zscore_random_sem = std(zscore_random, 0, 2, 'omitnan')/sqrt(nevent);
    zscore_mean = mean(zscore, 2, 'omitnan');
    zscore_sem = std(zscore, 0, 2, 'omitnan')/sqrt(nevent);
    
    hold on;
    patch([t; t(end:-1:1)], [zscore_random_mean+zscore_random_sem; zscore_random_mean(end:-1:1)-zscore_random_sem(end:-1:1)], [0.8 0.8 0.8], 'EdgeColor', [0.8 0.8 0.8]);
    plot(t, zscore_random_mean, '-', 'Color', [0 0 0]);
    patch([t; t(end:-1:1)], [zscore_mean+zscore_sem; zscore_mean(end:-1:1)-zscore_sem(end:-1:1)], [0 1 0], 'EdgeColor', [0 1 0]);
    plot(t, zscore_mean, '-', 'Color', [0 0.5 0]);
    yl = [min([zscore_random_mean-zscore_random_sem; zscore_mean-zscore_sem;]) max([zscore_random_mean+zscore_random_sem; zscore_mean+zscore_sem;])];
    line([0 0], yl, 'Color', [0 0 0], 'LineStyle', '--', 'LineWidth', 1);
    line([t(1) t(end)], [0 0], 'Color', [0 0 0], 'LineStyle', '-', 'LineWidth', 1);
    xlim([t(1) t(end)]);
    ylim(yl);
    set(gca, 'FontSize', 12);
    xlabel('Time (s)');
    ylabel('Z score');
    title(title_text);
end