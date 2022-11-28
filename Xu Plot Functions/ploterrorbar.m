function ploterrorbar(x, y, pcolor)
x = rmmissing(squeeze(x));
meanx = mean(x);
errorx = std(x)/sqrt(numel(x));
if y ~= 0
    y = rmmissing(squeeze(y));
    meany = mean(y);
    errory = std(y)/sqrt(numel(y));
    errorbar(meanx, meany, errory, errory, errorx, errorx, '.', 'Color', pcolor, 'CapSize', 15, 'MarkerSize', 12);
else
    errorbar(meanx, 0, errorx, 'horizontal', '.', 'Color', pcolor, 'CapSize', 15, 'MarkerSize', 12);
end