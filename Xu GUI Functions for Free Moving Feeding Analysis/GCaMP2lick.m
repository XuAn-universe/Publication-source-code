function GCaMP2lick(app, Exp_Path)
%*---------------------------------------------------------------------*
% Huang Lab
% Duke University
% Author : Xu An, Jan 2022
% xu.an@duke.edu
% Version: 1.0
%*---------------------------------------------------------------------*
value = app.TrialsListBox.Value;
value = sort(value);
trials = numel(value);

pre = app.presEditField.Value;
post = app.postsEditField.Value;

zscore_RetrievalStart_all = [];
zscore_RetrievalStart_random_all = [];
ID_Sit = false(1, trials);

% get photometry data
try
    fpdata_all = load([Exp_Path(1:end-7) '\FPData.mat']);
    nchannel = size(fpdata_all.zsignal_all, 1);
    for j = 1:nchannel
        lgdtext{j} = ['Channel ' num2str(j)];
    end
catch
    errordlg('Fiber photometry data is missing!', 'Error');
    return;
end

for i = 1:trials
    try
        temp = load([Exp_Path '\LabelledEvents' num2str(value(i)) '.mat']);
        events = temp.LabelledEvents;
    catch
        errordlg(['Trial ' num2str(value(i)) ' is not labelled'], 'Error');
        return;
    end
    
    fpdata = fpdata_all.zsignal_all(:, value(i));
    fpdata_zsignal = [];
    for j = 1:nchannel
        fpdata_zsignal(:, j) = fpdata{j}(:, 2);
    end
    fpdata_t = fpdata{1}(:, 1);
    SampleRate = mean(diff(fpdata_t));
    [t_aligned, zscore_RetrievalStart, zscore_RetrievalStart_random] = AlignSignal2Event(fpdata_t, fpdata_zsignal, events.RetrievalStart, pre, post);
    zscore_RetrievalStart_all = [zscore_RetrievalStart_all zscore_RetrievalStart];
    zscore_RetrievalStart_random_all = [zscore_RetrievalStart_random_all zscore_RetrievalStart_random];
    if ~isempty(events.SitStart)
        ID_Sit(i) = true;
    end
end

for i = 1:nchannel
    figure;
    subplot(1, 3, 1);
    plot_gcamp_align2event(t_aligned, zscore_RetrievalStart_random_all(:, :, i), zscore_RetrievalStart_all(:, :, i), ['aligned to RetrievalStart (all): ' lgdtext{i}]);
    set(gca, 'ButtonDownFcn', @extract_figure);
    subplot(1, 3, 2);
    plot_gcamp_align2event(t_aligned, zscore_RetrievalStart_random_all(:, ~ID_Sit, i), zscore_RetrievalStart_all(:, ~ID_Sit, i), ['aligned to RetrievalStart (no sit): ' lgdtext{i}]);
    set(gca, 'ButtonDownFcn', @extract_figure);
    subplot(1, 3, 3);
    plot_gcamp_align2event(t_aligned, zscore_RetrievalStart_random_all(:, ID_Sit, i), zscore_RetrievalStart_all(:, ID_Sit, i), ['aligned to RetrievalStart (with sit): ' lgdtext{i}]);
    set(gca, 'ButtonDownFcn', @extract_figure);
end

for i = 1:nchannel
    figure;
    heatmap4gcamp(t_aligned, zscore_RetrievalStart_all(:, :, i), ['aligned to RetrievalStart (all): ' lgdtext{i}]);
    if any(ID_Sit)
        figure;
        heatmap4gcamp(t_aligned, zscore_RetrievalStart_all(:, ID_Sit, i), ['aligned to RetrievalStart (with sit): ' lgdtext{i}]);
    end
    if any(~ID_Sit)
        figure;
        heatmap4gcamp(t_aligned, zscore_RetrievalStart_all(:, ~ID_Sit, i), ['aligned to RetrievalStart (no sit): ' lgdtext{i}]);
    end
end

colors = [0.3010 0.7450 0.9330; 0.8500 0.3250 0.0980];
hp = zeros(1, nchannel);
figure;
subplot(1, 3, 1);
for i = 1:nchannel
    hp(i) = plot_tj_MeanSEM(t_aligned, zscore_RetrievalStart_all(:, :, i), colors(i, :), colors(i, :), 'Time (s)', 'Z-score', 'aligned to RetrievalStart (all)');
end
yl = ylim;
plot([0 0], yl, '--k', 'LineWidth', 1);
legend(hp, lgdtext, 'Location', 'NorthEast');
legend('boxoff');

if any(~ID_Sit)
    subplot(1, 3, 2);
    for i = 1:nchannel
        hp(i) = plot_tj_MeanSEM(t_aligned, zscore_RetrievalStart_all(:, ~ID_Sit, i), colors(i, :), colors(i, :), 'Time (s)', 'Z-score', 'aligned to RetrievalStart (no sit)');
    end
    yl = ylim;
    plot([0 0], yl, '--k', 'LineWidth', 1);
    legend(hp, lgdtext, 'Location', 'NorthEast');
    legend('boxoff');
end

if any(ID_Sit)
    subplot(1, 3, 3);
    for i = 1:nchannel
        hp(i) = plot_tj_MeanSEM(t_aligned, zscore_RetrievalStart_all(:, ID_Sit, i), colors(i, :), colors(i, :), 'Time (s)', 'Z-score', 'aligned to RetrievalStart (with sit)');
    end
    yl = ylim;
    plot([0 0], yl, '--k', 'LineWidth', 1);
    legend(hp, lgdtext, 'Location', 'NorthEast');
    legend('boxoff');
end

result.zscore_lickstart = zscore_RetrievalStart_all(:, ~ID_Sit, :);
result.t = t_aligned;
assignin('base', 'result', result);