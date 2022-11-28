function plot_tj_individuals(t, data, individual_color, mean_color, xlabel_text, ylabel_text, title_text, filename)
if nargin < 8
    filename = [];
    if nargin < 7
        title_text = [];
    end
end
data_mean = mean(data, 2, 'omitnan');
hold on;
for i = 1:size(data, 2)
    if ~isempty(filename)
        plot(t, data(:, i), '-', 'Color', individual_color, 'ButtonDownFcn', @displayfilename, 'UserData', filename{i});
    else
        plot(t, data(:, i), '-', 'Color', individual_color, 'ButtonDownFcn', @displayfilename, 'UserData', []);
    end
end
plot(t, data_mean, '-', 'Color', mean_color, 'LineWidth', 2);
[xl(1), xl(2)] = bounds(t);
xlim(xl);
% [yl(1), yl(2)] = bounds(data(:));
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
set(gca, 'FontSize', 12, 'ButtonDownFcn', @extract_figure);

function displayfilename(src, eventdata)
if numel(src.UserData) > 4
    htext = text(eventdata.IntersectionPoint(1), eventdata.IntersectionPoint(2), src.UserData(1:end-4));
else
    htext = text(eventdata.IntersectionPoint(1), eventdata.IntersectionPoint(2), src.UserData);
end
pause(2);
try
    delete(htext)
    clear htext;
end