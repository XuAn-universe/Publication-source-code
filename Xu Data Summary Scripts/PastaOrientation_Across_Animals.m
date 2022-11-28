%*---------------------------------------------------------------------*
% Huang Lab
% Duke University
% Author : Xu An, Nov 2021
% xu.an@duke.edu
% Version: 1.0
%*---------------------------------------------------------------------*
%% pasta orientation analysis across mice (entire inhibition window)
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
FrameRate = 120;
trange = FrameRate*0.5;
stp = 2;
edges = 0:stp:90;
diffedges = -90:stp:90;

show_histogram = 1;

orientation_all_NI = [];
orientation_all_I = [];
orientationbite_all_NI = [];
orientationbite_all_I = [];
orientation_avg_NI = nan(1, N);
orientation_avg_I = nan(1, N);
orientation_var_NI = nan(1, N);
orientation_var_I = nan(1, N);
orientationbite_avg_NI = nan(1, N);
orientationbite_avg_I = nan(1, N);
p_all_NI = nan(N, numel(diffedges)-1);
p_all_I = nan(size(p_all_NI));
pbite_all_NI = nan(size(p_all_NI));
pbite_all_I = nan(size(p_all_NI));
Edist_biteNI = nan(N, 2);
correlationdist_biteNI = nan(N, 2);
spearmandist_biteNI = nan(N, 2);
Hdist_biteNI = nan(N, 2); % Hellinger distance
KLdivergence_biteNI = nan(N, 2); % Kullback–Leibler divergence
Edist_biteI = nan(N, 2);
correlationdist_biteI = nan(N, 2);
spearmandist_biteI = nan(N, 2);
Hdist_biteI = nan(N, 2); % Hellinger distance
KLdivergence_biteI = nan(N, 2); % Kullback–Leibler divergence

for i = 1:N
    temp = load([pathname filename{i}]);
    result = temp.result;
    temp = result.orientationxylight_all_NI{1};
    temp = rmmissing(temp);
    orientation_NI = temp;
    temp = result.orientationxylight_all_I{1};
    temp = rmmissing(temp);
    orientation_I = temp;
    
    temp = result.orientationxybite_NI{1};
    temp = temp(trange+1, :);
    temp = rmmissing(temp);
    orientationbite_NI = temp;
    temp = result.orientationxybite_I{1};
    if ~isempty(temp)
        temp = temp(trange+1, :);
        temp = rmmissing(temp);
        orientationbite_I = temp;
    else
        orientationbite_I = [];
    end
    
    if show_histogram
        figure;
        subplot(2, 1, 1);
        histogram(orientation_NI, edges, 'FaceColor', [0 0 0], 'FaceAlpha', 1, 'Normalization', 'probability');
        hold on;
        histogram(orientation_I, edges, 'FaceColor', 'none', 'EdgeColor', [0 1 0], 'FaceAlpha', 1, 'Normalization', 'probability');
        set(gca, 'TickDir', 'out', 'FontSize', 12, 'ButtonDownFcn', @extract_figure);
        box off;
        xlabel(['Orientation (' char(176) ')']);
        ylabel('probability');
        title('Orientation');
        
        subplot(2, 1, 2);
        histogram(orientationbite_NI, edges, 'FaceColor', [0 0 0], 'FaceAlpha', 1, 'Normalization', 'probability');
        if ~isempty(orientationbite_I)
            hold on;
            histogram(orientationbite_I, edges, 'FaceColor', 'none', 'EdgeColor', [0 1 0], 'FaceAlpha', 1, 'Normalization', 'probability');
        end
        set(gca, 'FontSize', 12, 'ButtonDownFcn', @extract_figure);
        box off;
        xlabel(['Orientation (' char(176) ')']);
        ylabel('probability');
        title('Orientation@bite');
    end
    
    [p_NI, ~] = histcounts(orientation_NI, edges, 'Normalization', 'probability');
    [p_I, ~] = histcounts(orientation_I, edges, 'Normalization', 'probability');
    [pbite_NI, ~] = histcounts(orientationbite_NI, edges, 'Normalization', 'probability');
    if ~isempty(orientationbite_I)
        [pbite_I, ~] = histcounts(orientationbite_I, edges, 'Normalization', 'probability');
    else
        pbite_I = [];
    end
    
    Edist_biteNI(i, 1) = pdist([p_NI;pbite_NI]);
    Edist_biteNI(i, 2) = pdist([p_I;pbite_NI]);
    
    correlationdist_biteNI(i, 1) = pdist([p_NI;pbite_NI], 'correlation');
    correlationdist_biteNI(i, 2) = pdist([p_I;pbite_NI], 'correlation');
    
    spearmandist_biteNI(i, 1) = pdist([p_NI;pbite_NI], 'spearman');
    spearmandist_biteNI(i, 2) = pdist([p_I;pbite_NI], 'spearman');
    
    Hdist_biteNI(i, 1) = 1/sqrt(2)*sqrt(sum((sqrt(p_NI)-sqrt(pbite_NI)).^2));
    Hdist_biteNI(i, 2) = 1/sqrt(2)*sqrt(sum((sqrt(p_I)-sqrt(pbite_NI)).^2));
    
    % the calculation below most likely is not meaningful as pbite_NI has a lot of zeros
    pbite_NIeps = pbite_NI;
    pbite_NIeps(pbite_NIeps == 0) = eps;
    temp = p_NI.*log2(p_NI./(pbite_NIeps));
    temp(isnan(temp)) = 0;
    KLdivergence_biteNI(i, 1) = sum(temp);
    temp = p_I.*log2(p_I./(pbite_NIeps));
    temp(isnan(temp)) = 0;
    KLdivergence_biteNI(i, 2) = sum(temp);
    
    if ~isempty(pbite_I)
        Edist_biteI(i, 1) = pdist([p_NI;pbite_I]);
        Edist_biteI(i, 2) = pdist([p_I;pbite_I]);
        
        correlationdist_biteI(i, 1) = pdist([p_NI;pbite_I], 'correlation');
        correlationdist_biteI(i, 2) = pdist([p_I;pbite_I], 'correlation');
        
        spearmandist_biteI(i, 1) = pdist([p_NI;pbite_I], 'spearman');
        spearmandist_biteI(i, 2) = pdist([p_I;pbite_I], 'spearman');
        
        Hdist_biteI(i, 1) = 1/sqrt(2)*sqrt(sum((sqrt(p_NI)-sqrt(pbite_I)).^2));
        Hdist_biteI(i, 2) = 1/sqrt(2)*sqrt(sum((sqrt(p_I)-sqrt(pbite_I)).^2));
        
        % the calculation below most likely is not meaningful as pbite_NI has a lot of zeros
        pbite_Ieps = pbite_I;
        pbite_Ieps(pbite_Ieps == 0) = eps;
        temp = p_NI.*log2(p_NI./(pbite_Ieps));
        temp(isnan(temp)) = 0;
        KLdivergence_biteI(i, 1) = sum(temp);
        temp = p_I.*log2(p_I./(pbite_Ieps));
        temp(isnan(temp)) = 0;
        KLdivergence_biteI(i, 2) = sum(temp);
    end
    
    orientation_all_NI = [orientation_all_NI; orientation_NI-mean(orientationbite_NI)];
    orientation_all_I = [orientation_all_I; orientation_I-mean(orientationbite_NI)];
    orientationbite_all_NI = [orientationbite_all_NI orientationbite_NI-mean(orientationbite_NI)];
    if ~isempty(orientationbite_I)
        orientationbite_all_I = [orientationbite_all_I orientationbite_I-mean(orientationbite_NI)];
    end
    orientation_avg_NI(i) = mean(orientation_NI);
    orientation_avg_I(i) = mean(orientation_I);
    orientationbite_avg_NI(i) = mean(orientationbite_NI);
    if ~isempty(orientationbite_I)
        orientationbite_avg_I(i) = mean(orientationbite_I);
    end
    
    orientation_var_NI(i) = var(orientation_NI);
    orientation_var_I(i) = var(orientation_I);
    
    [p_all_NI(i, :), ~] = histcounts(orientation_NI-mean(orientationbite_NI), diffedges, 'Normalization', 'probability');
    [p_all_I(i, :), ~] = histcounts(orientation_I-mean(orientationbite_NI), diffedges, 'Normalization', 'probability');
    [pbite_all_NI(i, :), ~] = histcounts(orientationbite_NI-mean(orientationbite_NI), diffedges, 'Normalization', 'probability');
    if ~isempty(orientationbite_I)
        [pbite_all_I(i, :), ~] = histcounts(orientationbite_I-mean(orientationbite_NI), diffedges, 'Normalization', 'probability');
    end
end
figure;
subplot(1, 5, 1);
paired_plot(Edist_biteNI(:, 1), Edist_biteNI(:, 2), 'Euclidean Distance', filename, 'bar');

subplot(1, 5, 2);
paired_plot(correlationdist_biteNI(:, 1), correlationdist_biteNI(:, 2), 'Correlation Distance', filename, 'bar');

subplot(1, 5, 3);
paired_plot(spearmandist_biteNI(:, 1), spearmandist_biteNI(:, 2), 'Spearman Distance', filename, 'bar');

subplot(1, 5, 4);
paired_plot(Hdist_biteNI(:, 1), Hdist_biteNI(:, 2), 'Hellinger Distance', filename, 'bar');

subplot(1, 5, 5);
paired_plot(KLdivergence_biteNI(:, 1), KLdivergence_biteNI(:, 2), 'Kullback–Leibler divergence', filename, 'bar');

figure;
subplot(1, 5, 1);
paired_plot(Edist_biteI(:, 1), Edist_biteI(:, 2), 'Euclidean Distance', filename, 'bar');

subplot(1, 5, 2);
paired_plot(correlationdist_biteI(:, 1), correlationdist_biteI(:, 2), 'Correlation Distance', filename, 'bar');

subplot(1, 5, 3);
paired_plot(spearmandist_biteI(:, 1), spearmandist_biteI(:, 2), 'Spearman Distance', filename, 'bar');

subplot(1, 5, 4);
paired_plot(Hdist_biteI(:, 1), Hdist_biteI(:, 2), 'Hellinger Distance', filename, 'bar');

subplot(1, 5, 5);
paired_plot(KLdivergence_biteI(:, 1), KLdivergence_biteI(:, 2), 'Kullback–Leibler divergence', filename, 'bar');

x = -90+stp/2:stp:90-stp/2;
figure;
plot_tj_MeanSEM(x', p_all_NI', [0 0 0], [0 0 0], [], [], []);
plot_tj_MeanSEM(x', p_all_I', [0 0 1], [0 0 1], ['Orientation (' char(176) ')'], 'Probability', 'Orientation');

compare2distributions(orientation_all_NI, orientation_all_I, ['Orientation (' char(176) ')'], 'Orientation');

figure;
plot_tj_MeanSEM(x', pbite_all_NI', [0 0 0], [0 0 0], [], [], []);
plot_tj_MeanSEM(x', pbite_all_I', [0 0 1], [0 0 1], ['Orientation (' char(176) ')'], 'Probability', 'Orientation@bite');

compare2distributions(orientationbite_all_NI, orientationbite_all_I, ['Orientation (' char(176) ')'], 'Orientation@bite');

figure;
plot(x, mean(p_all_NI), '-k');
hold on;
plot(x, mean(p_all_I), '-', 'Color', [0 0 0.5]);
plot(x, mean(pbite_all_NI), '-', 'Color', [0.8 0.8 0.8]);
plot(x, mean(pbite_all_I, 'omitnan'), '-b');
xlabel(['Orientation (' char(176) ')']);
ylabel('Probability');
set(gca, 'FontSize', 12);
box off;
title('Mean of individual mice');

[p_NI, ~] = histcounts(orientation_all_NI, diffedges, 'Normalization', 'probability');
[p_I, ~] = histcounts(orientation_all_I, diffedges, 'Normalization', 'probability');
[pbite_NI, ~] = histcounts(orientationbite_all_NI, diffedges, 'Normalization', 'probability');
[pbite_I, ~] = histcounts(orientationbite_all_I, diffedges, 'Normalization', 'probability');

figure;
plot(x, p_NI, '-k');
hold on;
plot(x, p_I, '-', 'Color', [0 0 0.5]);
plot(x, pbite_NI, '-', 'Color', [0.8 0.8 0.8]);
plot(x, pbite_I, '-b');
xlabel(['Orientation (' char(176) ')']);
ylabel('Probability');
set(gca, 'FontSize', 12);
box off;
title('Data of all mice combined');

figure;
paired_plot(orientation_avg_NI, orientation_avg_I, ['Orientation (' char(176) ')'], filename, 'bar');
figure;
paired_plot(orientationbite_avg_NI, orientationbite_avg_I, ['Orientation (' char(176) ')'], filename, 'bar');

figure;
paired_plot(orientation_var_NI, orientation_var_I, ['Orientation variance (' char(176) '^{2})'], filename, 'bar');

cd(curpwd);

%% pasta orientation at hand adjustment analysis across mice
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
FrameRate = 120;
trange = FrameRate*0.5;
stp = 2;
edges = 0:stp:90;
diffedges = -90:stp:90;

show_histogram = 1;

orientationadj_all_NI = [];
orientationadj_all_I = [];
orientationadj_avg_NI = nan(1, N);
orientationadj_avg_I = nan(1, N);
orientationadj_var_NI = nan(1, N);
orientationadj_var_I = nan(1, N);
padj_all_NI = nan(N, numel(diffedges)-1);
padj_all_I = nan(size(padj_all_NI));

for i = 1:N
    temp = load([pathname filename{i}]);
    result = temp.result;
    
    temp = result.orientationxyadj_NI{1};
    temp = temp(trange+1, :);
    temp = rmmissing(temp);
    orientationadj_NI = temp;
    temp = result.orientationxyadj_I{1};
    if ~isempty(temp)
        temp = temp(trange+1, :);
        temp = rmmissing(temp);
        orientationadj_I = temp;
    else
        orientationadj_I = [];
    end
    
    if show_histogram
        figure;
        histogram(orientationadj_NI, edges, 'FaceColor', [0 0 0], 'FaceAlpha', 1, 'Normalization', 'probability');
        if ~isempty(orientationadj_I)
            hold on;
            histogram(orientationadj_I, edges, 'FaceColor', 'none', 'EdgeColor', [0 1 0], 'FaceAlpha', 1, 'Normalization', 'probability');
        end
        set(gca, 'FontSize', 12, 'ButtonDownFcn', @extract_figure);
        box off;
        xlabel(['Orientation (' char(176) ')']);
        ylabel('probability');
        title('Orientation@adj');
    end

    orientationadj_all_NI = [orientationadj_all_NI orientationadj_NI-mean(orientationadj_NI)];
    if ~isempty(orientationadj_I)
        orientationadj_all_I = [orientationadj_all_I orientationadj_I-mean(orientationadj_NI)];
    end
    orientationadj_avg_NI(i) = mean(orientationadj_NI);
    orientationadj_var_NI(i) = var(orientationadj_NI);
    if ~isempty(orientationadj_I)
        orientationadj_avg_I(i) = mean(orientationadj_I);
        orientationadj_var_I(i) = var(orientationadj_I);
    end
    [padj_all_NI(i, :), ~] = histcounts(orientationadj_NI-mean(orientationadj_NI), diffedges, 'Normalization', 'probability');
    if ~isempty(orientationadj_I)
        [padj_all_I(i, :), ~] = histcounts(orientationadj_I-mean(orientationadj_NI), diffedges, 'Normalization', 'probability');
    end
end

x = -90+stp/2:stp:90-stp/2;

figure;
plot_tj_MeanSEM(x', padj_all_NI', [0 0 0], [0 0 0], [], [], []);
plot_tj_MeanSEM(x', padj_all_I', [0 0 1], [0 0 1], ['Orientation (' char(176) ')'], 'Probability', 'Orientation@adj');

compare2distributions(orientationadj_all_NI, orientationadj_all_I, ['Orientation (' char(176) ')'], 'Orientation@adj');

figure;
plot(x, mean(padj_all_NI), '-', 'Color', [0 0 0]);
hold on;
plot(x, mean(padj_all_I, 'omitnan'), '-', 'Color', [0.8 0.8 0.8]);
xlabel(['Orientation (' char(176) ')']);
ylabel('Probability');
set(gca, 'FontSize', 12);
box off;
title('Mean of individual mice');

[padj_NI, ~] = histcounts(orientationadj_all_NI, diffedges, 'Normalization', 'probability');
[padj_I, ~] = histcounts(orientationadj_all_I, diffedges, 'Normalization', 'probability');

figure;
plot(x, padj_NI, '-', 'Color', [0 0 0]);
hold on;
plot(x, padj_I, '-', 'Color', [0.8 0.8 0.8]);
xlabel(['Orientation (' char(176) ')']);
ylabel('Probability');
set(gca, 'FontSize', 12);
box off;
title('Data of all mice combined');

figure;
paired_plot(orientationadj_avg_NI, orientationadj_avg_I, ['Orientation (' char(176) ')'], filename, 'bar');

figure;
paired_plot(orientationadj_var_NI, orientationadj_var_I, ['Orientation variance (' char(176) '^{2})'], filename, 'bar');

cd(curpwd);