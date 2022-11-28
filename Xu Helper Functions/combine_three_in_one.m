%% Emx
h = findobj(gca, 'Type', 'line');
difE = nan(1, numel(h));
for i = 1:numel(h)
    difE(i) = diff(h(i).YData);
end

%% Fezf2
h = findobj(gca, 'Type', 'line');
difF = nan(1, numel(h));
for i = 1:numel(h)
    difF(i) = diff(h(i).YData);
end

%% PlexinD1
h = findobj(gca, 'Type', 'line');
difP = nan(1, numel(h));
for i = 1:numel(h)
    difP(i) = diff(h(i).YData);
end

%% bar plot
figure;
hold on;
bar(1, mean(difE), 0.8, 'FaceColor', [0 0 0]);
bar(2, mean(difF), 0.8, 'FaceColor', [0 1 0]);
bar(3, mean(difP), 0.8, 'FaceColor', [1 1 0]);
errorbar([1 2 3], [mean(difE) mean(difF) mean(difP)],...
    [std(difE)/sqrt(numel(difE)) std(difF)/sqrt(numel(difF)) std(difP)/sqrt(numel(difP))],...
    'k', 'LineStyle', 'none', 'CapSize', 15);
for i = 1:numel(difE)
    plot(1, difE(i), 'o', 'MarkerFaceColor', [1 1 1], 'MarkerEdgeColor', [0 0 0], 'MarkerSize', 10);
end
for i = 1:numel(difF)
    plot(2, difF(i), 'o', 'MarkerFaceColor', [1 1 1], 'MarkerEdgeColor', [0 0 0], 'MarkerSize', 10);
end
for i = 1:numel(difP)
    plot(3, difP(i), 'o', 'MarkerFaceColor', [1 1 1], 'MarkerEdgeColor', [0 0 0], 'MarkerSize', 10);
end
set(gca, 'XTick', [1 2 3], 'XTickLabel', {'Emx1', 'Fezf2', 'PlxnD1'}, 'XLim', [0.4 3.6], 'TickLength', [0 0], 'FontSize', 12);
ylabel('Test');

%% box plot
figure;
hold on;
group = [zeros(1, numel(difE)) ones(1, numel(difF)) ones(1, numel(difP))*2];
boxplot([difE difF difP], group, 'PlotStyle', 'traditional', 'Labels', {'Emx1', 'Fezf2', 'PlxnD1'}, 'Color', [0 0 0], 'LabelOrientation', 'horizontal', 'widths', 0.8);
for i = 1:numel(difE)
    plot(1, difE(i), 'o', 'MarkerFaceColor', [1 1 1], 'MarkerEdgeColor', [0 0 0], 'MarkerSize', 10);
end
for i = 1:numel(difF)
    plot(2, difF(i), 'o', 'MarkerFaceColor', [1 1 1], 'MarkerEdgeColor', [0 0 0], 'MarkerSize', 10);
end
for i = 1:numel(difP)
    plot(3, difP(i), 'o', 'MarkerFaceColor', [1 1 1], 'MarkerEdgeColor', [0 0 0], 'MarkerSize', 10);
end
set(gca, 'XLim', [0.4 3.6], 'FontSize', 12);
ylabel('Test');
box off;

%% box plot with 2 y axes
figure;
hold on;
group = [zeros(1, numel(difE)) ones(1, numel(difF)) ones(1, numel(difP))*2];
yyaxis left;
boxplot(difE, group(1:numel(difE)), 'PlotStyle', 'traditional', 'Labels', {'Emx1'}, 'Color', [0 0 0], 'LabelOrientation', 'horizontal', 'widths', 0.8);
for i = 1:numel(difE)
    plot(1, difE(i), 'o', 'MarkerFaceColor', [1 1 1], 'MarkerEdgeColor', [0 0 0], 'MarkerSize', 10);
end
ylim([8 28]);
ylabel('Test');
yyaxis right;
boxplot([difE difF difP], group, 'PlotStyle', 'traditional', 'Labels', {'Emx1', 'Fezf2', 'PlxnD1'}, 'Color', [1 0 0], 'LabelOrientation', 'horizontal', 'widths', 0.8);
for i = 1:numel(difF)
    plot(2, difF(i), 'o', 'MarkerFaceColor', [1 1 1], 'MarkerEdgeColor', [1 0 0], 'MarkerSize', 10);
end
for i = 1:numel(difP)
    plot(3, difP(i), 'o', 'MarkerFaceColor', [1 1 1], 'MarkerEdgeColor', [1 0 0], 'MarkerSize', 10);
end
ylim([-1 10]);
set(gca, 'XLim', [0.4 3.6], 'FontSize', 12);
box off;

%% 2 in 1 box plot
figure;
hold on;
group = [zeros(1, numel(difF)) ones(1, numel(difP))];
boxplot([difF difP], group, 'PlotStyle', 'traditional', 'Labels', {'Fezf2', 'PlxnD1'}, 'Color', [0 0 0], 'LabelOrientation', 'horizontal', 'widths', 0.8);
for i = 1:numel(difF)
    plot(1, difF(i), 'o', 'MarkerFaceColor', [1 1 1], 'MarkerEdgeColor', [0 0 0], 'MarkerSize', 10);
end
for i = 1:numel(difP)
    plot(2, difP(i), 'o', 'MarkerFaceColor', [1 1 1], 'MarkerEdgeColor', [0 0 0], 'MarkerSize', 10);
end
set(gca, 'XLim', [0.4 2.6], 'FontSize', 12);
ylabel('Test');
box off;