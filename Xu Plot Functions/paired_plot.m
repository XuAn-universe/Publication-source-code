function paired_plot(data_control, data_inhibition, ylabel_text, filename, style, x)
if nargin < 6
    x = [1 2];
end
if nargin < 5
    style = 'bar';
end
IDnan = isnan(data_control) | isnan(data_inhibition);
data_control(IDnan) = [];
data_inhibition(IDnan) = [];
if ~isempty(filename)
    filename(IDnan) = [];
end
if isempty(data_control) || isempty(data_inhibition)
    return;
end
hold on;
switch style
    case 'bar'
        linecolor = [0 0 0];
        bar(x(1), mean(data_control), 0.8, 'FaceColor', [0.5 0.5 0.5]);
        bar(x(2), mean(data_inhibition), 0.8, 'FaceColor', [0 1 0]);
        errorbar(x, [mean(data_control) mean(data_inhibition)],...
            [std(data_control)/sqrt(numel(data_control)) std(data_inhibition)/sqrt(numel(data_inhibition))],...
            'k', 'LineStyle', 'none', 'CapSize', 15);
        for i = 1:numel(data_control)
            if ~isempty(filename)
                plot(x, [data_control(i) data_inhibition(i)], '-o', 'Color', linecolor, 'MarkerFaceColor', [1 1 1], 'MarkerEdgeColor', linecolor, 'MarkerSize', 10, 'ButtonDownFcn', @displayfilename, 'UserData', filename{i});
            else
                plot(x, [data_control(i) data_inhibition(i)], '-o', 'Color', linecolor, 'MarkerFaceColor', [1 1 1], 'MarkerEdgeColor', linecolor, 'MarkerSize', 10, 'ButtonDownFcn', @displayfilename, 'UserData', []);
            end
        end
    case 'dot'
        linecolor = [0.8 0.8 0.8];
        for i = 1:numel(data_control)
            if ~isempty(filename)
                plot(x, [data_control(i) data_inhibition(i)], '-o', 'Color', linecolor, 'MarkerFaceColor', linecolor, 'MarkerEdgeColor', 'none', 'ButtonDownFcn', @displayfilename, 'UserData', filename{i});
            else
                plot(x, [data_control(i) data_inhibition(i)], '-o', 'Color', linecolor, 'MarkerFaceColor', linecolor, 'MarkerEdgeColor', 'none', 'ButtonDownFcn', @displayfilename, 'UserData', []);
            end
        end
        plot(x, [mean(data_control) mean(data_inhibition)], '-ok', 'MarkerFaceColor', [0 0 0], 'MarkerEdgeColor', 'none', 'LineWidth', 2);
        errorbar(x, [mean(data_control) mean(data_inhibition)],...
            [std(data_control)/sqrt(numel(data_control)) std(data_inhibition)/sqrt(numel(data_inhibition))],...
            'k', 'LineStyle', 'none', 'CapSize', 15);
end
set(gca, 'XTick', x, 'XTickLabel', {'Control', 'Inhibition'}, 'XLim', [x(1)-0.6 x(2)+0.6], 'FontSize', 12, 'ButtonDownFcn', @extract_figure);
ylabel(ylabel_text);
yl = ylim;
% ylim([0 yl(2)]);
hold off;

disp(['N = ' num2str(numel(data_control))]);
p = signrank(data_control, data_inhibition, 'tail', 'both');
disp(['Signrank: p = ' num2str(p)]);
[~, p] = ttest(data_control, data_inhibition, 'tail', 'both');
disp(['Paired Ttest: p = ' num2str(p)]);
try
    [H, ~, ~] = swtest(data_control-data_inhibition);
    if H
        disp('Data didn''t pass normality test');
    else
        disp('Data passed normality test');
    end
catch
    disp('Couldn''t perform swtest')
end

function displayfilename(src, eventdata)
htext = text(eventdata.IntersectionPoint(1), eventdata.IntersectionPoint(2), src.UserData(1:end-4));
pause(2);
try
    delete(htext)
    clear htext;
end