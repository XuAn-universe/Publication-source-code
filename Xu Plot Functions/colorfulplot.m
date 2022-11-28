function colorfulplot(x, y, xlabel_text, ylabel_text, title_text, cmap_text)
nline = numel(x)-1;
eval(['cmap = ' cmap_text '(' num2str(nline) ');']);
for i = 1:nline
    if ~any(isnan([x(i); x(i+1); y(i); y(i+1)]))
        line([x(i); x(i+1)], [y(i); y(i+1)], 'LineStyle', '-', 'LineWidth', 1, 'Color', cmap(i, :));
    end
end
xlabel(xlabel_text);
ylabel(ylabel_text);
title(title_text);
box off;
set(gca, 'TickDir', 'out', 'FontSize', 12, 'ButtonDownFcn', @extract_figure);