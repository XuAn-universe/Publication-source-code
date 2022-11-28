function nanID = heatmap4gcamp(t, signal, title_text)
nanID = [];
for i = 1:size(signal, 2)
    if any(isnan(signal(:, i)))
        nanID = [nanID i];
    end
end
signal(:, nanID) = [];
[tmin, tminID] = min(abs(t));
imshow(signal', [mean2(signal)-3*std2(signal) mean2(signal)+3*std2(signal)]);
colormap parula;
colorbar;
axis on;
hold on;
plot([tminID tminID], [1-0.5 size(signal, 2)+0.5], '--k', 'LineWidth', 1);
if size(signal, 2) == 1
    set(gca, 'XTick', [1 tminID numel(t)], 'XTickLabel', {num2str(t(1)), num2str(tmin), num2str(t(end))}, 'YTick', 1, 'YDir', 'normal');
else
    set(gca, 'XTick', [1 tminID numel(t)], 'XTickLabel', {num2str(t(1)), num2str(tmin), num2str(t(end))}, 'YTick', [1 size(signal, 2)], 'YDir', 'normal');
end
axis normal;
box off;
title(title_text);