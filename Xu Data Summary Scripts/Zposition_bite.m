%*---------------------------------------------------------------------*
% Huang Lab
% Duke University
% Author : Xu An, Oct 2021
% xu.an@duke.edu
% Version: 1.0
%*---------------------------------------------------------------------*
%%
curpwd = pwd;
try
   cd(pathname); 
end
[filename_saline, pathname] = uigetfile('*.mat', 'Pick all of the saline data set', 'MultiSelect', 'on');
if isequal(filename_saline, 0)
    cd(curpwd);
    return;
end
N = numel(filename_saline);
zbite_avg_saline = nan(1, N);
zbite_all_saline = cell(1, N);
zbite_saline = [];
for i = 1:N
    temp = load([pathname filename_saline{i}]);
    result = temp.result;
    temp = result.zbite_NI_1p{1};
    temp = rmmissing(temp);
    temp(temp > 20) = [];
    zbite_avg_saline(i) = mean(temp);
    zbite_all_saline{i} = temp;
    zbite_saline = [zbite_saline temp];
end

cd(pathname); 
[filename_muscimol, pathname] = uigetfile('*.mat', 'Pick all of the muscimol data set', 'MultiSelect', 'on');
if isequal(filename_muscimol, 0)
    cd(curpwd);
    return;
end
N = numel(filename_muscimol);
zbite_avg_muscimol = nan(1, N);
zbite_muscimol = [];
for i = 1:N
    temp = load([pathname filename_muscimol{i}]);
    result = temp.result;
    temp = result.zbite_NI_1p{1};
    temp = rmmissing(temp);
    temp(temp > 20) = [];
    zbite_avg_muscimol(i) = mean(temp);
    zbite_all_muscimol{i} = temp;
    zbite_muscimol = [zbite_muscimol temp];
end

cd(pathname); 
[filename_ref, pathname] = uigetfile('*.mat', 'Pick all of the reference data set', 'MultiSelect', 'on');
if isequal(filename_ref, 0)
    cd(curpwd);
    return;
end
N = numel(filename_ref);
zbite_avg_ref = nan(1, N);
zbite_ref = [];
for i = 1:N
    temp = load([pathname filename_ref{i}]);
    result = temp.result;
    temp = result.zbite_NI_1p{1};
    temp = rmmissing(temp);
    temp(temp > 20) = [];
    zbite_avg_ref(i) = mean(temp);
    zbite_all_ref{i} = temp;
    zbite_ref = [zbite_ref temp];
end

[~, edges] = histcounts([zbite_saline(:); zbite_muscimol(:)]);
figure;
histogram(zbite_saline, edges, 'FaceColor', [0 0 0], 'FaceAlpha', 1, 'Normalization', 'probability');
hold on;
histogram(zbite_muscimol, edges, 'FaceColor', 'none', 'EdgeColor', [0 1 0], 'FaceAlpha', 1, 'Normalization', 'probability');
yl = ylim;
line([mean(zbite_ref) mean(zbite_ref)], yl, 'Color', [1 0 0], 'LineWidth', 1, 'LineStyle', '--');
xl = xlim;
set(gca, 'TickDir', 'out', 'FontSize', 12, 'ButtonDownFcn', @extract_figure);
box off;
xlabel('Z (mm)');
ylabel('probability');
title('Data of all mice combined');

figure;
h1 = cdfplot(zbite_saline);
set(h1, 'Color', [0 0 0], 'LineWidth', 1.5);
hold on;
h2 = cdfplot(zbite_muscimol);
set(h2, 'Color', [0 1 0], 'LineWidth', 1.5);
xlim(xl);
set(gca, 'GridLineStyle', 'none', 'TickDir', 'out', 'FontSize', 12, 'ButtonDownFcn', @extract_figure);
box off;
xlabel('Z (mm)');
ylabel('Cumulative Probability');
[~, p] = kstest2(zbite_saline, zbite_muscimol);
disp(['KStest: p = ' num2str(p)]);
title('Data of all mice combined');

[~, edges] = histcounts([zbite_saline(:); zbite_muscimol(:); zbite_ref(:)]);
[p_saline, ~] = histcounts(zbite_saline, edges, 'Normalization', 'probability');
[p_muscimol, ~] = histcounts(zbite_muscimol, edges, 'Normalization', 'probability');
[p_ref, ~] = histcounts(zbite_ref, edges, 'Normalization', 'probability');
x = edges(1)+mean(diff(edges))/2:mean(diff(edges)):edges(end)-mean(diff(edges))/2;
figure;
plot(x, p_saline, '-k');
hold on;
plot(x, p_muscimol, '-g')
plot(x, p_ref, '-r')
set(gca, 'TickDir', 'out', 'FontSize', 12, 'ButtonDownFcn', @extract_figure);
box off;
xlabel('Z (mm)');
ylabel('probability');
title('Data of all mice combined');

figure;
h1 = cdfplot(zbite_saline);
set(h1, 'Color', [0 0 0], 'LineWidth', 1.5);
hold on;
h2 = cdfplot(zbite_muscimol);
set(h2, 'Color', [0 1 0], 'LineWidth', 1.5);
h3 = cdfplot(zbite_ref);
set(h3, 'Color', [1 0 0], 'LineWidth', 1.5);
xlim([min(edges) max(edges)]);
set(gca, 'GridLineStyle', 'none', 'TickDir', 'out', 'FontSize', 12, 'ButtonDownFcn', @extract_figure);
box off;
xlabel('Z (mm)');
ylabel('Cumulative Probability');
title('Data of all mice combined');
[~, p] = kstest2(zbite_saline, zbite_muscimol);
disp(['KStest: p = ' num2str(p)]);

p_all_saline = nan(N, numel(edges)-1);
p_all_muscimol = nan(size(p_all_saline));
p_all_ref = nan(size(p_all_saline));
for i = 1:N
    figure;
    histogram(zbite_all_saline{i}, edges, 'FaceColor', [0 0 0], 'FaceAlpha', 1, 'Normalization', 'probability');
    hold on;
    histogram(zbite_all_muscimol{i}, edges, 'FaceColor', 'none', 'EdgeColor', [0 1 0], 'FaceAlpha', 1, 'Normalization', 'probability');
    set(gca, 'TickDir', 'out', 'FontSize', 12, 'ButtonDownFcn', @extract_figure);
    box off;
    xlabel('Z (mm)');
    ylabel('probability');
    
    [p_all_saline(i, :), ~] = histcounts(zbite_all_saline{i}, edges, 'Normalization', 'probability');
    [p_all_muscimol(i, :), ~] = histcounts(zbite_all_muscimol{i}, edges, 'Normalization', 'probability');
    [p_all_ref(i, :), ~] = histcounts(zbite_all_ref{i}, edges, 'Normalization', 'probability');
end
figure;
plot_tj_MeanSEM(x', p_all_saline', [0 0 0], [0 0 0], [], [], []);
plot_tj_MeanSEM(x', p_all_muscimol', [0 0 1], [0 0 1], [], [], []);
plot_tj_MeanSEM(x', p_all_ref', [0 1 0], [0 1 0], 'Z (mm)', 'Probability', []);

% resampling = 10000;
% figure;
% permutation_test(zbite_all_saline, zbite_all_muscimol, 'both', resampling, 'Z (mm)');
% figure;
% bootstrapping_test(zbite_all_saline, zbite_all_muscimol, 'both', resampling, 'Z (mm)');

figure;
paired_plot(zbite_avg_saline, zbite_avg_muscimol, 'Z (mm)', filename_saline, 'dot');

violin_boxplot(zbite_saline, zbite_muscimol, 'Z (mm)');

cd(curpwd);