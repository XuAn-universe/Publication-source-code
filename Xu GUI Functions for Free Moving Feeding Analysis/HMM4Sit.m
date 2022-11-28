function HMM4Sit(app, Exp_Path, FrameRate, nmedian)
%*---------------------------------------------------------------------*
% Huang Lab
% Duke University
% Author : Xu An, Jan 2022
% xu.an@duke.edu
% Version: 1.0
%*---------------------------------------------------------------------*
exctime = 0.5;
value = app.TrialsListBox.Value;
value = sort(value);
trials = numel(value);

FrameRate = round(FrameRate);

sitstart_labelled = cell(1, trials);
zdata = cell(1, trials);
zdata_all = [];
time4z = cell(1, trials);

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
    [~, ~, ~, ~, ~, ~, ~, znose, ~, ~, ~, ~] = trajectory_postprocessing(3, Exp_Path, value(i), label_table, nmedian, FrameRate);
    [~, ~, ~, ~, ~, ~, ~, zpawl, ~, ~, ~, ~] = trajectory_postprocessing(10, Exp_Path, value(i), label_table, nmedian, FrameRate);
    [~, ~, ~, ~, ~, ~, ~, zpawr, ~, ~, ~, ~] = trajectory_postprocessing(16, Exp_Path, value(i), label_table, nmedian, FrameRate);
    
    t = (1:size(znose, 1))'/FrameRate;
    
    bite_timestamps = Bite_events(value(i)).time_bites;
    znose(t < exctime | t >= bite_timestamps(1), :) = [];
    znose = fillmissing(znose(:, 1), 'spline');
    zpawl(t < exctime | t >= bite_timestamps(1), :) = [];
    zpawl = fillmissing(zpawl(:, 1), 'spline');
    zpawr(t < exctime | t >= bite_timestamps(1), :) = [];
    zpawr = fillmissing(zpawr(:, 1), 'spline');
    t(t < exctime | t >= bite_timestamps(1)) = [];
    zdata_all = [zdata_all [znose'; zpawl'; zpawr']];
    zdata{i} = [znose'; zpawl'; zpawr'];
    time4z{i} = t;
    
    try
        temp = load([Exp_Path '\LabelledEvents' num2str(value(i)) '.mat']);
        events = temp.LabelledEvents;
        temp = events.SitStart;
        temp(temp > bite_timestamps(1)) = [];
        temp = temp(end); % only consider the sit right before 1st bite
        sitstart_labelled{i} = temp;
    catch
        sitstart_labelled{i} = [];
    end
end
% figure;
% subplot(1, 3, 1);
% histogram(zdata_all(1, :), 'FaceColor', [0 0 0], 'FaceAlpha', 1, 'Normalization', 'probability');
% set(gca, 'TickDir', 'out', 'FontSize', 12, 'ButtonDownFcn', @extract_figure);
% box off;
% xlabel('Z nose (mm)');
% ylabel('probability');
% subplot(1, 3, 2);
% histogram(zdata_all(2, :), 'FaceColor', [0 0 0], 'FaceAlpha', 1, 'Normalization', 'probability');
% set(gca, 'TickDir', 'out', 'FontSize', 12, 'ButtonDownFcn', @extract_figure);
% box off;
% xlabel('Z left paw (mm)');
% ylabel('probability');
% subplot(1, 3, 3);
% histogram(zdata_all(3, :), 'FaceColor', [0 0 0], 'FaceAlpha', 1, 'Normalization', 'probability');
% set(gca, 'TickDir', 'out', 'FontSize', 12, 'ButtonDownFcn', @extract_figure);
% box off;
% xlabel('Z right paw (mm)');
% ylabel('probability');

ninitialization = 10;
O = size(zdata_all, 1); % number of signals
Q = 2; % number of states
M = 1; % number of Gaussians (per state per signal)

cov_type = 'full'; % covariance of Gaussians

[mu0, Sigma0] = mixgauss_init(Q*M, zdata_all, cov_type);
mu0 = reshape(mu0, [O Q M]);
Sigma0 = reshape(Sigma0, [O O Q M]);

evalmatrix = nan(trials, ninitialization);
for j = 1:ninitialization
    if j == 1
        % initial guess of parameters
        prior0 = [1 0]; % always start with locomotion, assign locomotion state 1
        transmat0 = [0.99 0.01; 0.01 0.99]; % based on real data
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
        mhmm_em(zdata, prior0, transmat0, mu0, Sigma0, mixmat0, 'max_iter', 1000);
    
    if M ~= 1
        state1mean = sum(squeeze(mu1(:, 1, :)).*repmat(mixmat1(1, :), O, 1), 'all');
        state2mean = sum(squeeze(mu1(:, 2, :)).*repmat(mixmat1(2, :), O, 1), 'all');
    else
        state1mean = sum(mu1(:, 1, 1));
        state2mean = sum(mu1(:, 2, 1));
    end
    if state1mean <= state2mean
        swaplabels = 0;
    else
        swaplabels = 1;
    end

    for i = 1:trials
        B = mixgauss_prob(zdata{i}, mu1, Sigma1, mixmat1);
        statepath = viterbi_path(prior1, transmat1, B);
        if swaplabels
            temp = statepath;
            statepath(temp == 1) = 2;
            statepath(temp == 2) = 1;
        end
        gttime = sitstart_labelled{i}; %ground truth timestamp of sit start
        t = time4z{i};
        figure;
        plotstate(t, statepath, [min(zdata{i}, [], 'all') max(zdata{i}, [], 'all')], [1 0 0; 0 1 0]);
        hold on;
        for k = 1:O
            plot(t, zdata{i}(k, :), '-', 'Color', ones(1, 3)*(k-1)/O);
        end
        if ~isempty(gttime)
            plot([gttime gttime], [min(zdata{i}, [], 'all') max(zdata{i}, [], 'all')], '--k');
            sitID = find(diff(statepath) == 1);
            if ~isempty(sitID)
                predsittime = t(sitID(end)+1);
                evalmatrix(i, j) = predsittime-gttime;
            end
        end
    end
end

figure;
plot(LL, '-ok');
set(gca, 'TickDir', 'out', 'FontSize', 12, 'ButtonDownFcn', @extract_figure);
box off;
xlabel('Iteration');
ylabel('Log-likelihood');