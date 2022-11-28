function violin_boxplot(data_control, data_inhibition, ylabel_text)
data_control = data_control(:)';
data_inhibition = data_inhibition(:)';
figure;
subplot(1, 2, 1);
violin({data_control data_inhibition}, 'facecolor', [0.8 0.8 0.8; 0 1 0],...
    'facealpha', 1, 'mc', [], 'medc', []);
yl = ylim;
hold on;
group = [zeros(1, numel(data_control)) ones(1, numel(data_inhibition))];
boxplot([data_control data_inhibition], group, 'PlotStyle', 'compact', 'colors', [0 0 0], 'Symbol', 'o',...
    'Labels', {'Control', 'Inhibition'}, 'LabelOrientation', 'horizontal', 'widths', 0.8);
ylabel(ylabel_text);
box off;
ylim(yl);
hold off;
set(gca, 'FontSize', 12, 'ButtonDownFcn', @extract_figure);

subplot(1, 2, 2);
boxplot([data_control data_inhibition], group, 'PlotStyle', 'traditional', 'Labels', {'Control', 'Inhibition'}, 'LabelOrientation', 'horizontal', 'widths', 0.8);
ylabel(ylabel_text);
box off;
set(gca, 'FontSize', 12, 'ButtonDownFcn', @extract_figure);

p = ranksum(data_control, data_inhibition, 'tail', 'both');
disp(['Ranksum: p = ' num2str(p) ' for ' ylabel_text]);
[~, p] = ttest2(data_control, data_inhibition, 'tail', 'both', 'Vartype', 'equal');
disp(['Ttest (equal unknown variance): p = ' num2str(p) ' for ' ylabel_text]);
[~, p] = ttest2(data_control, data_inhibition, 'tail', 'both', 'Vartype', 'unequal');
disp(['Ttest (unequal unknown variance): p = ' num2str(p) ' for ' ylabel_text]);
[H, ~, ~] = swtest(data_control);
if H
    disp('Control didn''t pass normality test');
else
    disp('Control passed normality test');
end
[H, ~, ~] = swtest(data_inhibition);
if H
    disp('Inhibition didn''t pass normality test');
else
    disp('Inhibition passed normality test');
end