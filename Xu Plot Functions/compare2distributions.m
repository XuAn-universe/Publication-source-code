function compare2distributions(NI, I, xlabel_text, title_text)
[~, edges] = histcounts([NI(:); I(:)]);
[xs, xl] = bounds([NI(:); I(:)]);
figure;
subplot(3, 1, 1);
histogram(NI, edges, 'FaceColor', [0 0 0], 'FaceAlpha', 1, 'Normalization', 'probability');
hold on;
histogram(I, edges, 'FaceColor', 'none', 'EdgeColor', [0 1 0], 'FaceAlpha', 1, 'Normalization', 'probability');
set(gca, 'TickDir', 'out', 'FontSize', 12, 'ButtonDownFcn', @extract_figure);
box off;
xlabel(xlabel_text);
ylabel('probability');
title(title_text);

try
    logNI = log(NI);
    logI = log(I);
    [~, logedges] = histcounts([logNI(:); logI(:)]);
    subplot(3, 1, 2);
    histogram(logNI, logedges, 'FaceColor', [0 0 0], 'FaceAlpha', 1, 'Normalization', 'probability');
    hold on;
    histogram(logI, logedges, 'FaceColor', 'none', 'EdgeColor', [0 1 0], 'FaceAlpha', 1, 'Normalization', 'probability');
    set(gca, 'TickDir', 'out', 'FontSize', 12, 'ButtonDownFcn', @extract_figure);
    box off;
    xlabel(['log ' xlabel_text]);
    ylabel('probability');
    title(title_text);
end

subplot(3, 1, 3);
hold on;
h1 = cdfplot(NI);
set(h1, 'Color', [0 0 0], 'LineWidth', 1.5);
h2 = cdfplot(I);
set(h2, 'Color', [0 1 0], 'LineWidth', 1.5);
xlim([xs, xl]);
set(gca, 'GridLineStyle', 'none', 'TickDir', 'out', 'FontSize', 12, 'ButtonDownFcn', @extract_figure);
xlabel(xlabel_text);
ylabel('Cumulative Probability');
legend({'Control', 'Inhibition'}, 'FontSize', 12);
legend('boxoff');
title(title_text);

[~, p] = kstest2(NI, I, 'tail', 'unequal');
disp(['KStest: p = ' num2str(p) ' for ' xlabel_text]);