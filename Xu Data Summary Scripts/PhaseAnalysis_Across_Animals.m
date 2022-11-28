%*---------------------------------------------------------------------*
% Huang Lab
% Duke University
% Author : Xu An, July 2022
% xu.an@duke.edu
% Version: 1.0
%*---------------------------------------------------------------------*
%%
N = 2;
result.phase_NI = [];
result.phase_I = [];
for i = 1:N
    eval(['result.phase_NI = [result.phase_NI result' num2str(i) 's.phase_NI];']);
    eval(['result.phase_I = [result.phase_I result' num2str(i) 's.phase_I];']);
end

%%
N = 2;
result.phase_NI = [];
result.phase_I = [];
for i = 1:N
    eval(['result.phase_NI = [result.phase_NI result' num2str(i) 'g.phase_NI];']);
    eval(['result.phase_I = [result.phase_I result' num2str(i) 'g.phase_I];']);
end

%% bite phase analysis across mice
curpwd = pwd;
try
   cd(pathname); 
end
[filename, pathname] = uigetfile('*.mat', 'Pick all of the data set', 'MultiSelect', 'on');
if isequal(filename, 0)
    cd(curpwd);
    return;
end

N = numel(filename);

clc;
stp = 18;
edges = 0:stp:360;
check_individual = 1;

phase_all_NI = [];
phase_all_I = [];
phase_count_NI = nan(N, numel(edges)-1);
phase_count_I = nan(N, numel(edges)-1);
phase_avg_NI = nan(1, N);
phase_avg_I = nan(1, N);
r_avg_NI = nan(1, N);
r_avg_I = nan(1, N);

for i = 1:N
    temp = load([pathname filename{i}]);
    result = temp.result;
    phase_all_NI = [phase_all_NI result.phase_NI];
    phase_all_I = [phase_all_I result.phase_I];
    phase_count_NI(i, :) = histcounts(result.phase_NI, edges, 'Normalization', 'probability');
    phase_count_I(i, :) = histcounts(result.phase_I, edges, 'Normalization', 'probability');
    [mu, ul, ll] = circ_mean(circ_ang2rad(result.phase_NI'));
    phase_avg_NI(i) = mu;
    [mu, ul, ll] = circ_mean(circ_ang2rad(result.phase_I'));
    phase_avg_I(i) = mu;
    r_avg_NI(i) = circ_r(circ_ang2rad(result.phase_NI'));
    r_avg_I(i) = circ_r(circ_ang2rad(result.phase_I'));
    
    if check_individual
        figure('Name', filename{i});
        subplot(1, 2 ,1);
        histogram(result.phase_NI, edges, 'FaceColor', [0 0 0], 'EdgeColor', [0 0 0], 'FaceAlpha', 1, 'Normalization', 'probability');
        hold on;
        histogram(result.phase_I, edges, 'FaceColor', 'none', 'EdgeColor', [0 1 1], 'FaceAlpha', 1, 'Normalization', 'probability');
        plot(edges, max([phase_count_NI(i, :) phase_count_I(i, :)])/2*cosd(edges)+max([phase_count_NI(i, :) phase_count_I(i, :)])/2, '-r');
        xlim([0 360]);
        set(gca, 'XTick', 0:90:360, 'TickDir', 'out', 'FontSize', 12, 'ButtonDownFcn', @extract_figure);
        box off;
        xlabel(['Phase (' char(176) ')']);
        ylabel('Probability');
        
        subplot(1, 2 ,2);
        polarhistogram(circ_ang2rad(result.phase_NI), edges/180*pi, 'FaceColor', [0 0 0], 'EdgeColor', [0 0 0], 'FaceAlpha', 0.5, 'Normalization', 'probability');
        hold on;
        polarhistogram(circ_ang2rad(result.phase_I), edges/180*pi, 'FaceColor', [0 1 1], 'EdgeColor', [0 1 1], 'FaceAlpha', 0.5, 'Normalization', 'probability');
        set(gca, 'FontSize', 12, 'ButtonDownFcn', @extract_figure);
    end
end

figure;
paired_plot(phase_avg_NI/pi*180, phase_avg_I/pi*180, ['Phase (' char(176) ')'], filename, 'bar');
figure;
paired_plot(r_avg_NI, r_avg_I, 'Vector length', filename, 'bar');

figure;
histogram(phase_all_NI, edges, 'FaceColor', [0 0 0], 'EdgeColor', [0 0 0], 'FaceAlpha', 1, 'Normalization', 'probability');
hold on;
histogram(phase_all_I, edges, 'FaceColor', 'none', 'EdgeColor', [0 1 1], 'FaceAlpha', 1, 'Normalization', 'probability');
% plot(edges, max([phase_count_NI(i, :) phase_count_I(i, :)])/2*cosd(edges)+max([phase_count_NI(i, :) phase_count_I(i, :)])/2, '-r');
xlim([0 360]);
set(gca, 'XTick', 0:90:360, 'TickDir', 'out', 'FontSize', 12, 'ButtonDownFcn', @extract_figure);
box off;
xlabel(['Phase (' char(176) ')']);
ylabel('Probability');

[pval, table] = circ_wwtest(circ_ang2rad(phase_all_NI), circ_ang2rad(phase_all_I));
disp(['WW test: p = ' num2str(pval)]);
[pval, med, P] = circ_cmtest(circ_ang2rad(phase_all_NI), circ_ang2rad(phase_all_I));
disp(['CM test: p = ' num2str(pval)]);
[pval, k, K] = circ_kuipertest(circ_ang2rad(phase_all_NI), circ_ang2rad(phase_all_I));
disp(['Kuiper test: p < ' num2str(pval)]);

figure;
hp(1) = plot_tj_MeanSEM(stp/2:stp:360-stp/2, phase_count_NI', [0 0 0], [0 0 0], ['Phase (' char(176) ')'], 'Probability', '');
hp(2) = plot_tj_MeanSEM(stp/2:stp:360-stp/2, phase_count_I', [0 0 1], [0 0 1], ['Phase (' char(176) ')'], 'Probability', '');
xlim([0 360]);
set(gca, 'XTick', 0:90:360, 'TickDir', 'out', 'FontSize', 12, 'ButtonDownFcn', @extract_figure);
legend(hp, {'Ctr', 'Inh'}, 'FontSize', 12);
legend('boxoff');

cd(curpwd);

%%
x = [1 2];
data_control = [phase_avg_NIs phase_avg_NIg]/pi*180;
data_inhibition = [phase_avg_Is phase_avg_Ig]/pi*180;
figure;
hold on;
bar(x(1), mean(data_control), 0.8, 'FaceColor', [0.5 0.5 0.5]);
bar(x(2), mean(data_inhibition), 0.8, 'FaceColor', [0 1 0]);
errorbar(x, [mean(data_control) mean(data_inhibition)],...
    [std(data_control)/sqrt(numel(data_control)) std(data_inhibition)/sqrt(numel(data_inhibition))],...
    'k', 'LineStyle', 'none', 'CapSize', 15);
linecolor = [1 0 0];
for i = 1:numel(phase_avg_NIs)
    plot(x, [phase_avg_NIs(i) phase_avg_Is(i)]/pi*180, '-o', 'Color', linecolor, 'MarkerFaceColor', [1 1 1], 'MarkerEdgeColor', linecolor, 'MarkerSize', 10, 'ButtonDownFcn', @displayfilename, 'UserData', []);
end
linecolor = [0 0 1];
for i = 1:numel(phase_avg_NIs)
    plot(x, [phase_avg_NIg(i) phase_avg_Ig(i)]/pi*180, '-o', 'Color', linecolor, 'MarkerFaceColor', [1 1 1], 'MarkerEdgeColor', linecolor, 'MarkerSize', 10, 'ButtonDownFcn', @displayfilename, 'UserData', []);
end
set(gca, 'XTick', x, 'XTickLabel', {'Control', 'Inhibition'}, 'XLim', [x(1)-0.6 x(2)+0.6], 'FontSize', 12, 'ButtonDownFcn', @extract_figure);
ylabel(['Phase (' char(176) ')']);
figure;
paired_plot([phase_avg_NIs phase_avg_NIg]/pi*180, [phase_avg_Is phase_avg_Ig]/pi*180, ['Phase (' char(176) ')'], '', 'bar');

data_control = [r_avg_NIs r_avg_NIg];
data_inhibition = [r_avg_Is r_avg_Ig];
figure;
hold on;
bar(x(1), mean(data_control), 0.8, 'FaceColor', [0.5 0.5 0.5]);
bar(x(2), mean(data_inhibition), 0.8, 'FaceColor', [0 1 0]);
errorbar(x, [mean(data_control) mean(data_inhibition)],...
    [std(data_control)/sqrt(numel(data_control)) std(data_inhibition)/sqrt(numel(data_inhibition))],...
    'k', 'LineStyle', 'none', 'CapSize', 15);
linecolor = [1 0 0];
for i = 1:numel(r_avg_NIs)
    plot(x, [r_avg_NIs(i) r_avg_Is(i)], '-o', 'Color', linecolor, 'MarkerFaceColor', [1 1 1], 'MarkerEdgeColor', linecolor, 'MarkerSize', 10, 'ButtonDownFcn', @displayfilename, 'UserData', []);
end
linecolor = [0 0 1];
for i = 1:numel(r_avg_NIs)
    plot(x, [r_avg_NIg(i) r_avg_Ig(i)], '-o', 'Color', linecolor, 'MarkerFaceColor', [1 1 1], 'MarkerEdgeColor', linecolor, 'MarkerSize', 10, 'ButtonDownFcn', @displayfilename, 'UserData', []);
end
set(gca, 'XTick', x, 'XTickLabel', {'Control', 'Inhibition'}, 'XLim', [x(1)-0.6 x(2)+0.6], 'FontSize', 12, 'ButtonDownFcn', @extract_figure);
ylabel('Vector length');
figure;
paired_plot([r_avg_NIs r_avg_NIg], [r_avg_Is r_avg_Ig], 'Vector length', '', 'bar');

%% bite phase analysis across mice for control condition
curpwd = pwd;
try
   cd(pathname); 
end
[filename, pathname] = uigetfile('*.mat', 'Pick all of the data set', 'MultiSelect', 'on');
if isequal(filename, 0)
    cd(curpwd);
    return;
end

N = numel(filename);

clc;
stp = 18;
edges = 0:stp:360;
check_individual = 1;

phase_all_NI = [];
phase_count_NI = nan(N, numel(edges)-1);
phase_avg_NI = nan(1, N);
r_avg_NI = nan(1, N);

for i = 1:N
    temp = load([pathname filename{i}]);
    result = temp.result;
    phase_all_NI = [phase_all_NI result.phase_NI];
    phase_count_NI(i, :) = histcounts(result.phase_NI, edges, 'Normalization', 'probability');
    [mu, ul, ll] = circ_mean(circ_ang2rad(result.phase_NI'));
    phase_avg_NI(i) = mu;
    r_avg_NI(i) = circ_r(circ_ang2rad(result.phase_NI'));
    
    if check_individual
        figure('Name', filename{i});
        subplot(1, 2 ,1);
        histogram(result.phase_NI, edges, 'FaceColor', [0 0 0], 'EdgeColor', [0 0 0], 'FaceAlpha', 1, 'Normalization', 'probability');
        hold on;
        plot(edges, max(phase_count_NI(i, :))/2*cosd(edges)+max(phase_count_NI(i, :))/2, '-r');
        xlim([0 360]);
        set(gca, 'XTick', 0:90:360, 'TickDir', 'out', 'FontSize', 12, 'ButtonDownFcn', @extract_figure);
        box off;
        xlabel(['Phase (' char(176) ')']);
        ylabel('Probability');
        
        subplot(1, 2 ,2);
        polarhistogram(circ_ang2rad(result.phase_NI), edges/180*pi, 'FaceColor', [0 0 0], 'EdgeColor', [0 0 0], 'FaceAlpha', 0.5, 'Normalization', 'probability');
        set(gca, 'FontSize', 12, 'ButtonDownFcn', @extract_figure);
    end
end

figure;
histogram(phase_all_NI, edges, 'FaceColor', [0 0 0], 'EdgeColor', [0 0 0], 'FaceAlpha', 1, 'Normalization', 'probability');
xlim([0 360]);
set(gca, 'XTick', 0:90:360, 'TickDir', 'out', 'FontSize', 12, 'ButtonDownFcn', @extract_figure);
box off;
xlabel(['Phase (' char(176) ')']);
ylabel('Probability');

figure;
hp(1) = plot_tj_MeanSEM(stp/2:stp:360-stp/2, phase_count_NI', [0 0 0], [0 0 0], ['Phase (' char(176) ')'], 'Probability', '');
xlim([0 360]);
set(gca, 'XTick', 0:90:360, 'TickDir', 'out', 'FontSize', 12, 'ButtonDownFcn', @extract_figure);
legend(hp, {'Ctr'}, 'FontSize', 12);
legend('boxoff');

cd(curpwd);

%%
figure;
paired_plot(phase_avg_NIs/pi*180, phase_avg_NIg/pi*180, ['Phase (' char(176) ')'], '', 'bar');
figure;
paired_plot(r_avg_NIs, r_avg_NIg, 'Vector length', '', 'bar');