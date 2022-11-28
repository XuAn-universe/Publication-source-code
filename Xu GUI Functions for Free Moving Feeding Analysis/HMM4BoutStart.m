function HMM4BoutStart(app, Exp_Path, FrameRate, nmedian)
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

FrameRate = round(FrameRate);
chewdurationmin = 0.1; % threshold for chewing duration

boutstart_labelled = cell(1, trials);
distancedata = cell(1, trials);
distancedata_all = [];
time4distance = cell(1, trials);

if app.TrackingDataCheckBox.Value
    label_table = table2array(app.UITable.Data);
else
    helpdlg('''Tracking Data'' needs to be checked');
    return;
end

try
    audiolocation = Exp_Path(1:end-7);
    temp = load([audiolocation '\Detected_Bite_Events.mat']);
    Bite_events = temp.Audio_analysis;
catch
    errordlg('Please detect the bites first.', 'Error');
    return;
end

for i = 1:trials
    try
        temp = load([Exp_Path '\LabelledEvents' num2str(value(i)) '.mat']);
        events = temp.LabelledEvents;
        boutstart_labelled{i} = events.BiteBoutStart;
    catch
        boutstart_labelled{i} = [];
    end
    
    [~, ~, ~, ~, ~, ~, ~, znose, ~, ~, ~, ~] = trajectory_postprocessing(3, Exp_Path, value(i), label_table, nmedian, FrameRate);
    [~, ~, ~, ~, ~, ~, ~, zpawl, ~, ~, ~, ~] = trajectory_postprocessing(10, Exp_Path, value(i), label_table, nmedian, FrameRate);
    [~, ~, ~, ~, ~, ~, ~, zpawr, ~, ~, ~, ~] = trajectory_postprocessing(16, Exp_Path, value(i), label_table, nmedian, FrameRate);
    
    t = (1:size(znose, 1))'/FrameRate;
    
    bite_timestamps = Bite_events(value(i)).time_bites;
%     bite_timestamps(1) = 10; %
%     bite_timestamps(end) = 65; %
    pawl2nose = zpawl(:, 1)-znose(:, 1);
    pawl2nose = fillmissing(pawl2nose, 'spline');
    pawl2nose(t < bite_timestamps(1) | t > bite_timestamps(end)) = [];
    pawl2nose = (pawl2nose-mean(pawl2nose)); % mean is a bit better than median
    pawr2nose = zpawr(:, 1)-znose(:, 1);
    pawr2nose = fillmissing(pawr2nose, 'spline');
    pawr2nose(t < bite_timestamps(1) | t > bite_timestamps(end)) = [];
    pawr2nose = (pawr2nose-mean(pawr2nose)); % mean is a bit better than median
    t(t < bite_timestamps(1) | t > bite_timestamps(end)) = [];
    distancedata_all = [distancedata_all [pawl2nose'; pawr2nose']];
    distancedata{i} = [pawl2nose'; pawr2nose'];
    time4distance{i} = t;
end
% figure;
% subplot(1, 2, 1);
% histogram(distancedata_all(1, :), 'FaceColor', [0 0 0], 'FaceAlpha', 1, 'Normalization', 'probability');
% set(gca, 'TickDir', 'out', 'FontSize', 12, 'ButtonDownFcn', @extract_figure);
% box off;
% xlabel('Left paw to nose distance (mm)');
% ylabel('probability');
% subplot(1, 2, 2);
% histogram(distancedata_all(2, :), 'FaceColor', [0 0 0], 'FaceAlpha', 1, 'Normalization', 'probability');
% set(gca, 'TickDir', 'out', 'FontSize', 12, 'ButtonDownFcn', @extract_figure);
% box off;
% xlabel('Right paw to nose distance (mm)');
% ylabel('probability');

ninitialization = 10;
O = size(distancedata_all, 1); % number of signals
Q = 2; % number of states
M = 1; % number of Gaussians (per state per signal)

cov_type = 'full'; % covariance of Gaussians

[mu0, Sigma0] = mixgauss_init(Q*M, distancedata_all, cov_type);
mu0 = reshape(mu0, [O Q M]);
Sigma0 = reshape(Sigma0, [O O Q M]);
figure;
for i = 1:ninitialization
    if i == 1
        % initial guess of parameters
        prior0 = [0 1]; % always start with handling, assign handling state 2
        transmat0 = [0.98 0.02; 0.02 0.98]; % based on real data
    else
        % initial guess of parameters
        prior0 = normalise(rand(Q,1));
        transmat0 = mk_stochastic(rand(Q,Q));
    end
    if M ~= 1
        mixmat0 = mk_stochastic(rand(Q,M));
    else
        mixmat0 = [];
    end
    
    [LL, prior1, transmat1, mu1, Sigma1, mixmat1] = ...
        mhmm_em(distancedata, prior0, transmat0, mu0, Sigma0, mixmat0, 'max_iter', 100);
    
    subplot(2, ceil(ninitialization/2), i);
    plot(LL, '-ok');
    set(gca, 'TickDir', 'out', 'FontSize', 12, 'ButtonDownFcn', @extract_figure);
    box off;
    xlabel('Iteration');
    ylabel('Log-likelihood');
    
    modelparameter(i).LL = LL(end);
    modelparameter(i).prior1 = prior1;
    modelparameter(i).transmat1 = transmat1;
    modelparameter(i).mu1 = mu1;
    modelparameter(i).Sigma1 = Sigma1;
    modelparameter(i).mixmat1 = mixmat1;
end
if app.EvaluateCheckBox.Value
    evalmatrix = nan(ninitialization, 3);
    terror_all = cell(1, ninitialization);
    terror_refined_all = cell(1, ninitialization);
    hf = figure;
    for j = 1:ninitialization
        prior1 = modelparameter(j).prior1;
        transmat1 = modelparameter(j).transmat1;
        mu1 = modelparameter(j).mu1;
        Sigma1 = modelparameter(j).Sigma1;
        mixmat1 = modelparameter(j).mixmat1;
        % check to see whether state labels need to be swapped
        if M ~= 1
            state1mean = mu1(1, 1, :).*mixmat1(1, :)+mu1(2, 1, :).*mixmat1(1, :);
            state2mean = mu1(1, 2, :).*mixmat1(2, :)+mu1(2, 2, :).*mixmat1(2, :);
        else
            state1mean = mu1(1, 1, :)+mu1(2, 1, :);
            state2mean = mu1(1, 2, :)+mu1(2, 2, :);
        end
        if state1mean <= state2mean
            swaplabels = 0;
        else
            swaplabels = 1;
        end
        nextra_all = 0;
        nhit_all = 0;
        nmiss_all = 0;
        terror_all{j} = [];
        terror_refined_all{j} = [];
        figure;
        for i = 1:trials
            B = mixgauss_prob(distancedata{i}, mu1, Sigma1, mixmat1);
            statepath = viterbi_path(prior1, transmat1, B);
            if swaplabels
                temp = statepath;
                statepath(temp == 1) = 2;
                statepath(temp == 2) = 1;
            end
            t = time4distance{i};
            raise_interval = find_raise_interval(statepath); % get the start and stop IDs of all chew-handle bouts, as well as IDs of handle bout start
            raise_interval(:, (t(raise_interval(3, :))-t(raise_interval(1, :))) < chewdurationmin) = [];
            boutstartHMM = compute_boutstartHMM(raise_interval, distancedata{i}); % refine the IDs of handle bout start based on raise height
            gttime = boutstart_labelled{i}; % ground truth timestamps
            if i <= 12
                subplot(3, 4, i);
                plotstate(t, statepath, [min(distancedata{i}, [], 'all') max(distancedata{i}, [], 'all')], [1 0 0; 0 1 0]);
                hold on;
                for k = 1:O
                    plot(t, distancedata{i}(k, :), '-', 'Color', ones(1, 3)*(k-1)/O);
                end
                plot([t(boutstartHMM)'; t(boutstartHMM)'], repmat([min(distancedata{i}, [], 'all'); max(distancedata{i}, [], 'all')], 1, numel(boutstartHMM)), '--b');
                %             plot(t(raise_interval(1, :)), ones(size(raise_interval, 2), 1)*max(distancedata{i}, [], 'all'), 'v', 'Color', [0.8 0.8 0.8]);
                %             plot(t(raise_interval(2, :)), ones(size(raise_interval, 2), 1)*max(distancedata{i}, [], 'all'), 'v', 'Color', [0 0 0]);
                if ~isempty(gttime)
                    plot([gttime'; gttime'], repmat([min(distancedata{i}, [], 'all'); max(distancedata{i}, [], 'all')], 1, numel(gttime)), '--k');
                end
                set(gca, 'TickDir', 'out', 'FontSize', 12, 'ButtonDownFcn', @extract_figure);
                box off;
                xlabel('Time (s)');
                ylabel('Distance (mm)');
            end
            if ~isempty(gttime)
                if ~isempty(raise_interval)
                    result = analyze_performance(gttime, raise_interval, t);
                    result_refined = analyze_performance(gttime, [raise_interval(1:2, :); boutstartHMM], t);
                    nextra_all = nextra_all+result.nextra;
                    nhit_all = nhit_all+result.nhit;
                    nmiss_all = nmiss_all+result.nmiss;
                    terror_all{j} = [terror_all{j}; result.terror];
                    terror_refined_all{j} = [terror_refined_all{j}; result_refined.terror];
                else
                    nmiss_all = nmiss_all+numel(gttime);
                end
            end
        end
        evalmatrix(j, :) = [nextra_all nhit_all nmiss_all];
        figure(hf);
        subplot(2, ceil(ninitialization/2), j);
        histogram(terror_all{j}, 'FaceColor', [0 0 0], 'EdgeColor', 'none', 'FaceAlpha', 1, 'Normalization', 'probability');
        hold on;
        histogram(terror_refined_all{j}, 'FaceColor', 'none', 'EdgeColor', [0 1 0], 'FaceAlpha', 1, 'Normalization', 'probability');
        yl = ylim;
        plot(ones(1, 2)*mean(abs(terror_all{j})), yl, '--', 'Color', [0 0 0]);
        plot(ones(1, 2)*mean(abs(terror_refined_all{j})), yl, '--', 'Color', [0 1 0]);
        set(gca, 'TickDir', 'out', 'FontSize', 12, 'ButtonDownFcn', @extract_figure);
        box off;
        xlabel('\Delta time (s)');
        ylabel('probability');
    end
    extrarate = evalmatrix(:, 1)./(evalmatrix(:, 1)+evalmatrix(:, 2))
    hitrate = evalmatrix(:, 2)./(evalmatrix(:, 3)+evalmatrix(:, 2))
    evalresult.evalmatrix = evalmatrix;
    evalresult.terror = terror_all;
    evalresult.terror_refined = terror_refined_all;
    assignin('base', 'result', evalresult);
else
    for i = 1:ninitialization
        if i == 1
            LLmax = modelparameter(1).LL;
            bestmodel = 1;
        else
            if modelparameter(i).LL > LLmax
                LLmax = modelparameter(i).LL;
                bestmodel = i;
            end
        end
    end
    prior1 = modelparameter(bestmodel).prior1;
    transmat1 = modelparameter(bestmodel).transmat1;
    mu1 = modelparameter(bestmodel).mu1;
    Sigma1 = modelparameter(bestmodel).Sigma1;
    mixmat1 = modelparameter(bestmodel).mixmat1;
    % check to see whether state labels need to be swapped
    if M ~= 1
        state1mean = mu1(1, 1, :).*mixmat1(1, :)+mu1(2, 1, :).*mixmat1(1, :);
        state2mean = mu1(1, 2, :).*mixmat1(2, :)+mu1(2, 2, :).*mixmat1(2, :);
    else
        state1mean = mu1(1, 1, :)+mu1(2, 1, :);
        state2mean = mu1(1, 2, :)+mu1(2, 2, :);
    end
    if state1mean <= state2mean
        swaplabels = 0;
    else
        swaplabels = 1;
    end
    figure;
    for i = 1:trials
        B = mixgauss_prob(distancedata{i}, mu1, Sigma1, mixmat1);
        statepath = viterbi_path(prior1, transmat1, B);
        if swaplabels
            temp = statepath;
            statepath(temp == 1) = 2;
            statepath(temp == 2) = 1;
        end
        t = time4distance{i};
        raise_interval = find_raise_interval(statepath); % get the start and stop IDs of all chew-handle bouts, as well as IDs of handle bout start
        raise_interval(:, (t(raise_interval(3, :))-t(raise_interval(1, :))) < chewdurationmin) = [];
        boutstartHMM = compute_boutstartHMM(raise_interval, distancedata{i}); % refine the IDs of handle bout start based on raise height
        if i <= 12
            if trials ~= 1
                subplot(3, 4, i);
            end
            plotstate(t, statepath, [min(distancedata{i}, [], 'all') max(distancedata{i}, [], 'all')], [1 0 0; 0 1 0]);
            hold on;
            for k = 1:O
                plot(t, distancedata{i}(k, :), '-', 'Color', ones(1, 3)*(k-1)/O);
            end
            plot([t(boutstartHMM)'; t(boutstartHMM)'], repmat([min(distancedata{i}, [], 'all'); max(distancedata{i}, [], 'all')], 1, numel(boutstartHMM)), '--b');
            set(gca, 'TickDir', 'out', 'FontSize', 12, 'ButtonDownFcn', @extract_figure);
            box off;
            xlabel('Time (s)');
            ylabel('Distance (mm)');
        end
        
        temp = load([Exp_Path '\LabelledEvents' num2str(value(i)) '.mat']);
        LabelledEvents = temp.LabelledEvents;
        LabelledEvents.ChewStartHMM = t(raise_interval(1, :));
        if raise_interval(1, 1) == 1
            LabelledEvents.ChewStartHMM(1) = [];
        end
%         LabelledEvents.BiteBoutStartHMM = t(raise_interval(3, :)); % uncomment it!
        save([Exp_Path '\LabelledEvents' num2str(value(i)) '.mat'], 'LabelledEvents');
    end
    if trials == 1
        assignin('base', 'LabelledEvents', LabelledEvents);
    end
end
msgbox('Done !');