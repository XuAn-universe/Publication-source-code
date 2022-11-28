function event4heatmap(t, tevent)
hold on;
for i = 1:numel(tevent)
    if ~isnan(tevent(i)) && tevent(i) <= max(t)
        [~, tminID] = min(abs(t-tevent(i)));
        plot([tminID tminID], [i-0.5 i+0.5], '-k', 'LineWidth', 1);
    end
end
set(gca, 'YDir', 'reverse', 'ButtonDownFcn', @extract_figure);