function histogram_multitrial(NI, bite_NI, I, bite_I, xlabel_text, title_text)
if ~isempty(I)
    [~, edges] = histcounts([NI(:); I(:)]);
else
    [~, edges] = histcounts(NI(:));
end

pd = fitdist(NI(:), 'normal');
y = pdf(pd, edges)*(edges(2)-edges(1));
lowbound = pd.mu-1*pd.sigma;
highbound = pd.mu+1*pd.sigma;
xmedian = median(NI, 'omitnan');

figure;
subplot(3, 1, 1);
histogram(NI, edges, 'FaceColor', [0 0 0], 'FaceAlpha', 1, 'Normalization', 'probability');
hold on;
if ~isempty(I)
    histogram(I, edges, 'FaceColor', 'none', 'EdgeColor', [0 1 0], 'FaceAlpha', 1, 'Normalization', 'probability');
end
plot(edges, y, '-r', 'LineWidth', 1);
yl = ylim;
line([lowbound lowbound], yl, 'Color', [1, 0, 0], 'LineWidth', 1, 'LineStyle', '--');
line([highbound highbound], yl, 'Color', [1, 0, 0], 'LineWidth', 1, 'LineStyle', '--');
line([xmedian xmedian], yl, 'Color', [0, 0, 1], 'LineWidth', 1, 'LineStyle', '--');
xl = xlim;
set(gca, 'TickDir', 'out', 'FontSize', 12, 'ButtonDownFcn', @extract_figure);
box off;
xlabel(xlabel_text);
ylabel('probability');
title(title_text);

subplot(3, 1, 3);
hold on;
h1 = cdfplot(NI);
set(h1, 'Color', [0 0 0], 'LineWidth', 1.5);

if ~isempty(bite_NI)
    h2 = cdfplot(bite_NI);
    set(h2, 'Color', [0.75 0.75 0.75], 'LineWidth', 1.5);
end

if ~isempty(I)
    h3 = cdfplot(I);
    set(h3, 'Color', [0 0.5 0], 'LineWidth', 1.5);
end

if ~isempty(bite_I)
    h4 = cdfplot(bite_I);
    set(h4, 'Color', [0 1 0], 'LineWidth', 1.5);
end

xlim(xl);
set(gca, 'GridLineStyle', ':', 'TickDir', 'out', 'FontSize', 12, 'ButtonDownFcn', @extract_figure);
xlabel(xlabel_text);
ylabel('Cumulative Probability');
if ~isempty(I) && ~isempty(bite_NI) && ~isempty(bite_I)
    legend({'Control', 'Control@bite', 'Inhibition', 'Inhibition@bite'}, 'FontSize', 12);
end
if isempty(bite_NI) && isempty(bite_I) && ~isempty(I)
    legend({'Control', 'Inhibition'}, 'FontSize', 12);
end
if isempty(I) && isempty(bite_I) && ~isempty(bite_NI)
    legend({'Control', 'Control@bite'}, 'FontSize', 12);
end
if isempty(bite_I) && ~isempty(I) && ~isempty(bite_NI)
    legend({'Control', 'Control@bite', 'Inhibition'}, 'FontSize', 12);
end
if isempty(bite_I) && isempty(I) && isempty(bite_NI)
    legend('Control', 'FontSize', 12);
end
legend('boxoff');
title(title_text);

if isempty(bite_NI) && isempty(bite_I)
    return;
end

if ~isempty(bite_I)
    [~, edges] = histcounts([bite_NI(:); bite_I(:)]);
else
    [~, edges] = histcounts(bite_NI(:));
end

pd = fitdist(bite_NI(:), 'normal');
y = pdf(pd, edges)*(edges(2)-edges(1));
lowbound = pd.mu-1*pd.sigma;
highbound = pd.mu+1*pd.sigma;
xmedian = median(bite_NI, 'omitnan');

subplot(3, 1, 2);
histogram(bite_NI, edges, 'FaceColor', [0 0 0], 'FaceAlpha', 1, 'Normalization', 'probability');
hold on;
if ~isempty(bite_I)
    histogram(bite_I, edges, 'FaceColor', 'none', 'EdgeColor', [0 1 0], 'FaceAlpha', 1, 'Normalization', 'probability');
end
plot(edges, y, '-r', 'LineWidth', 1);
yl = ylim;
line([lowbound lowbound], yl, 'Color', [1, 0, 0], 'LineWidth', 1, 'LineStyle', '--');
line([highbound highbound], yl, 'Color', [1, 0, 0], 'LineWidth', 1, 'LineStyle', '--');
line([xmedian xmedian], yl, 'Color', [0, 0, 1], 'LineWidth', 1, 'LineStyle', '--');
xlim(xl);
set(gca, 'TickDir', 'out', 'FontSize', 12, 'ButtonDownFcn', @extract_figure);
box off;
xlabel(xlabel_text);
ylabel('probability');
title([title_text '@bite']);