function bar_plot(data_control, data_inhibition, ylabel_text, x)
if nargin < 4
    x = [1 2];
end
hold on;
bar(x(1), mean(data_control), 0.8, 'FaceColor', [0.5 0.5 0.5]);
plot(x(1)+rand(1, numel(data_control))*0.8-0.4, data_control, 'o',...
    'MarkerFaceColor', [0 0 0], 'MarkerEdgeColor', 'none');
bar(x(2), mean(data_inhibition), 0.8, 'FaceColor', [0 1 0]);
plot(x(2)+rand(1, numel(data_inhibition))*0.8-0.4, data_inhibition, 'o',...
    'MarkerFaceColor', [0 0 0], 'MarkerEdgeColor', 'none');
errorbar(x, [mean(data_control) mean(data_inhibition)],...
    [std(data_control)/sqrt(numel(data_control)) std(data_inhibition)/sqrt(numel(data_inhibition))],...
    'k', 'LineStyle', 'none', 'CapSize', 15);
set(gca, 'XTick', x, 'XTickLabel', {'Control', 'Inhibition'}, 'XLim', [x(1)-0.6 x(2)+0.6], 'FontSize', 12);
ylabel(ylabel_text);
hold off;

p = ranksum(data_control, data_inhibition, 'tail', 'both');
disp(['Ranksum: p = ' num2str(p)]);
[~, p] = ttest2(data_control, data_inhibition, 'tail', 'both');
disp(['Ttest2: p = ' num2str(p)]);