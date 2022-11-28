function p = bootstrapping_test(x, y, tail, resampling, xlabel_text)
% x and y need to be paired
N = numel(x);
xmean_all = nan(1, N);
ymean_all = nan(1, N);
for i = 1:N
    xmean_all(i) = mean(x{i}(:));
    ymean_all(i) = mean(y{i}(:));
end
diffobs = mean(ymean_all)-mean(xmean_all);
diffdistr = nan(1, resampling);
for i = 1:resampling
    xmean_all = nan(1, N);
    ymean_all = nan(1, N);
    for j = 1:N
        nx = numel(x{j});
        ny = numel(y{j});
        xy = [x{j}(:); y{j}(:)];
        xmean_all(j) = mean(xy(randi(nx+ny, 1, nx)));
        ymean_all(j) = mean(xy(randi(nx+ny, 1, ny)));
    end
    diffdistr(i) = mean(ymean_all)-mean(xmean_all);
end
switch tail
    case 'left'
        p = (sum(diffdistr < diffobs)+1)/(resampling+1);
    case 'right'
        p = (sum(diffdistr > diffobs)+1)/(resampling+1);
    case 'both'
        p = (sum(abs(diffdistr) > abs(diffobs))+1)/(resampling+1);
end

histogram(diffdistr, 'FaceColor', [0 0 0], 'FaceAlpha', 1, 'Normalization', 'probability');
hold on;
yl = ylim;
line([diffobs diffobs], yl, 'Color', [1, 0, 0], 'LineWidth', 1, 'LineStyle', '--');
switch tail
    case 'both'
        line([-diffobs -diffobs], yl, 'Color', [1, 0, 0], 'LineWidth', 1, 'LineStyle', '--');
end
set(gca, 'TickLength', [0 0], 'FontSize', 12);
box off;
xlabel(xlabel_text);
ylabel('probability');
title(['Bootstrapping Test: p = ' num2str(p)]);
hold off;