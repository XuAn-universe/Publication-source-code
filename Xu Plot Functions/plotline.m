function plotline(x1, x2, y1, y2, lcolor)
x1mean = mean(x1, 'omitnan');
x2mean = mean(x2, 'omitnan');
y1mean = mean(y1, 'omitnan');
y2mean = mean(y2, 'omitnan');
plot([x1mean x2mean], [y1mean y2mean], '-', 'Color', lcolor);