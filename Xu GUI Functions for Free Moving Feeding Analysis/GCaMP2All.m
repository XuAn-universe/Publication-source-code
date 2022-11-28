function GCaMP2All(app, Exp_Path)
value = app.TrialsListBox.Value;
value = sort(value);
trials = numel(value);

pre = app.presEditField.Value;
post = app.postsEditField.Value;
gap = 0.1; % x axis gap between plots

zscore_stepcross_all = [];
zscore_stepcross_random_all = [];
zscore_retrievalstart_all = [];
zscore_retrievalstart_random_all = [];
zscore_withdraw_all = [];
zscore_withdraw_random_all = [];
zscore_handadj_all = [];
zscore_handadj_random_all = [];
zscore_bite_all = [];
zscore_bite_random_all = [];
zscore_chew_all = [];
zscore_chew_random_all = [];

retrievalstart2stepcross = [];

firstwithdraw2retrieval = [];
firstadj2retrieval = [];
firstbite2retrieval = [];

firstadj2withdraw = [];
firstbite2withdraw = [];
firstchew2withdraw = [];

firstbite2adj = [];

firstadj2bite = [];

load([Exp_Path '\Analysis_Session.mat'], 'Video_annotation');
    
Bite_events = [];
% get bite events
try
    audiolocation = Exp_Path(1:end-7);
    temp = load([audiolocation '\Detected_Bite_Events.mat']);
    Bite_events = temp.Audio_analysis;
end
    
% get photometry data
try
    fpdata_all = load([audiolocation '\FPData.mat']);
    nchannel = size(fpdata_all.zsignal_all, 1);
    for j = 1:nchannel
        lgdtext{j} = ['Channel ' num2str(j)];
    end
catch
    errordlg('Fiber photometry data is missing!', 'Error');
    return;
end

for i = 1:trials
    if ~Video_annotation(value(i)).Disgard
        try
            temp = load([Exp_Path '\LabelledEvents' num2str(value(i)) '.mat']);
            events = temp.LabelledEvents;
        catch
            continue;
        end
        
        bite_timestamps = [];
        if ~isempty(Bite_events)
            bite_timestamps = Bite_events(value(i)).time_bites;
            bite_timestamps = sort(bite_timestamps);
        end
        
        fpdata = fpdata_all.zsignal_all(:, value(i));
        fpdata_zsignal = [];
        fpdata_t = [];
        for j = 1:nchannel
            fpdata_zsignal(:, j) = fpdata{j}(:, 2);
        end
        fpdata_t = fpdata{1}(:, 1);
        SampleRate = mean(diff(fpdata_t));
   
        adjstart = [];
        if ~isempty(events.PawLAdjustmentStart) || ~isempty(events.PawRAdjustmentStart)
            adjstart = get_adjustment_start(events.PawRAdjustmentStart, events.PawRAdjustmentEnd, events.PawLAdjustmentStart, events.PawLAdjustmentEnd);
            adjstart = sort(adjstart);
            [t_aligned, zscore_handadj, zscore_handadj_random] = AlignSignal2Event(fpdata_t, fpdata_zsignal, adjstart, pre, post);
            zscore_handadj_all = [zscore_handadj_all zscore_handadj];
            zscore_handadj_random_all = [zscore_handadj_random_all zscore_handadj_random];
            
            for j = 1:numel(adjstart)
                biteafteradj = [];
                if ~isempty(bite_timestamps)
                    biteafteradj = bite_timestamps(bite_timestamps > adjstart(j));
                end
                if ~isempty(biteafteradj)
                    firstbite2adj = [firstbite2adj; biteafteradj(1)-adjstart(j)];
                else
                    firstbite2adj = [firstbite2adj; NaN];
                end
            end
        end
        
        if isfield(events, 'ChewStartHMM') && ~isempty(events.ChewStartHMM)
            events.ChewStartHMM = sort(events.ChewStartHMM);
            [t_aligned, zscore_chew, zscore_chew_random] = AlignSignal2Event(fpdata_t, fpdata_zsignal, events.ChewStartHMM, pre, post);
            zscore_chew_all = [zscore_chew_all zscore_chew];
            zscore_chew_random_all = [zscore_chew_random_all zscore_chew_random];
        end
        
        withdraw = [];
        if isfield(events, 'BiteBoutStartHMM') && ~isempty(events.BiteBoutStartHMM)
            withdraw = events.BiteBoutStartHMM;
        elseif ~isempty(events.BiteBoutStart)
            withdraw = events.BiteBoutStart;
        end
        if ~isempty(withdraw)
            withdraw = sort(withdraw);
            [t_aligned, zscore_withdraw, zscore_withdraw_random] = AlignSignal2Event(fpdata_t, fpdata_zsignal, withdraw, pre, post);
            zscore_withdraw_all = [zscore_withdraw_all zscore_withdraw];
            zscore_withdraw_random_all = [zscore_withdraw_random_all zscore_withdraw_random];
            
            for j = 1:numel(withdraw)
                adjafterwithdraw = [];
                if ~isempty(adjstart)
                    adjafterwithdraw = adjstart(adjstart > withdraw(j));
                end
                if ~isempty(adjafterwithdraw)
                    firstadj2withdraw = [firstadj2withdraw; adjafterwithdraw(1)-withdraw(j)];
                else
                    firstadj2withdraw = [firstadj2withdraw; NaN];
                end
                
                biteafterwithdraw = [];
                if ~isempty(bite_timestamps)
                    biteafterwithdraw = bite_timestamps(bite_timestamps > withdraw(j));
                end
                if ~isempty(biteafterwithdraw)
                    firstbite2withdraw = [firstbite2withdraw; biteafterwithdraw(1)-withdraw(j)];
                else
                    firstbite2withdraw = [firstbite2withdraw; NaN];
                end
                
                chewafterwithdraw = [];
                if isfield(events, 'ChewStartHMM') && ~isempty(events.ChewStartHMM)
                    chewafterwithdraw = events.ChewStartHMM(events.ChewStartHMM > withdraw(j));
                end
                if ~isempty(chewafterwithdraw)
                    firstchew2withdraw = [firstchew2withdraw; chewafterwithdraw(1)-withdraw(j)];
                else
                    firstchew2withdraw = [firstchew2withdraw; NaN];
                end
            end
        end
        
        if ~isempty(events.RetrievalStart)
            events.RetrievalStart = events.RetrievalStart(ismember(events.RetrievalStart, events.MouthRetrievalStart)); % only include successful retrieval using mouth
            events.RetrievalStart = sort(events.RetrievalStart);
            if ~isempty(events.RetrievalStart)
                [t_aligned, zscore_retrievalstart, zscore_retrievalstart_random] = AlignSignal2Event(fpdata_t, fpdata_zsignal, events.RetrievalStart, pre, post);
                zscore_retrievalstart_all = [zscore_retrievalstart_all zscore_retrievalstart];
                zscore_retrievalstart_random_all = [zscore_retrievalstart_random_all zscore_retrievalstart_random];
                
                [t_aligned, zscore_stepcross, zscore_stepcross_random] = AlignSignal2Event(fpdata_t, fpdata_zsignal, 0, pre, post);
                zscore_stepcross_all = [zscore_stepcross_all zscore_stepcross];
                zscore_stepcross_random_all = [zscore_stepcross_random_all zscore_stepcross_random];
                
                retrievalstart2stepcross = [retrievalstart2stepcross; min(events.RetrievalStart)];
            end
            
            for j = 1:numel(events.RetrievalStart)
                withdrawafterretrieval = [];
                if ~isempty(withdraw)
                    withdrawafterretrieval = withdraw(withdraw > events.RetrievalStart(j));
                end
                if ~isempty(withdrawafterretrieval)
                    firstwithdraw2retrieval = [firstwithdraw2retrieval; withdrawafterretrieval(1)-events.RetrievalStart(j)];
                else
                    firstwithdraw2retrieval = [firstwithdraw2retrieval; NaN];
                end
                
                adjafterretrieval = [];
                if ~isempty(adjstart)
                    adjafterretrieval = adjstart(adjstart > events.RetrievalStart(j));
                end
                if ~isempty(adjafterretrieval)
                    firstadj2retrieval = [firstadj2retrieval; adjafterretrieval(1)-events.RetrievalStart(j)];
                else
                    firstadj2retrieval = [firstadj2retrieval; NaN];
                end
                
                biteafterretrieval = [];
                if ~isempty(bite_timestamps)
                    biteafterretrieval = bite_timestamps(bite_timestamps > events.RetrievalStart(j));
                end
                if ~isempty(biteafterretrieval)
                    firstbite2retrieval = [firstbite2retrieval; biteafterretrieval(1)-events.RetrievalStart(j)];
                else
                    firstbite2retrieval = [firstbite2retrieval; NaN];
                end
            end
        end

        if ~isempty(bite_timestamps)
            [t_aligned, zscore_bite, zscore_bite_random] = AlignSignal2Event(fpdata_t, fpdata_zsignal, bite_timestamps, pre, post);
            zscore_bite_all = [zscore_bite_all zscore_bite];
            zscore_bite_random_all = [zscore_bite_random_all zscore_bite_random];
            
            for j = 1:numel(bite_timestamps)
                adjafterbite = [];
                if ~isempty(adjstart)
                    adjafterbite = adjstart(adjstart > bite_timestamps(j));
                end
                if ~isempty(adjafterbite)
                    firstadj2bite = [firstadj2bite; adjafterbite(1)-bite_timestamps(j)];
                else
                    firstadj2bite = [firstadj2bite; NaN];
                end
            end
        end
    end
end

for i = 1:nchannel
    if ~isempty(zscore_stepcross_all)
        figure('Name', ['aligned to Step Cross: ' lgdtext{i}]);
        subplot(1, 2, 1);
        heatmap4gcamp(t_aligned, zscore_stepcross_all(:, :, i), 'Sort by trial number');
        set(gca, 'ButtonDownFcn', @extract_figure);
        
        [retrievalstart2stepcross_sort, IDsort] = sort(retrievalstart2stepcross);
        IDsort(isnan(retrievalstart2stepcross_sort)) = [];
        subplot(1, 2, 2);
        nanID = heatmap4gcamp(t_aligned, zscore_stepcross_all(:, IDsort, i), 'Sort by Retrieval Start');
        event4heatmap(t_aligned, retrievalstart2stepcross_sort(setdiff(1:numel(IDsort), nanID)));
    end
    
    if ~isempty(zscore_retrievalstart_all)
        figure('Name', ['aligned to Retrieval Start with Mouth: ' lgdtext{i}]);
        subplot(2, 2, 1);
        heatmap4gcamp(t_aligned, zscore_retrievalstart_all(:, :, i), 'Sort by trial number');
        set(gca, 'ButtonDownFcn', @extract_figure);
        
        [firstwithdraw2retrieval_sort, IDsort] = sort(firstwithdraw2retrieval);
        IDsort(isnan(firstwithdraw2retrieval_sort)) = [];
        subplot(2, 2, 2);
        nanID = heatmap4gcamp(t_aligned, zscore_retrievalstart_all(:, IDsort, i), 'Sort by withdraw');
        event4heatmap(t_aligned, firstwithdraw2retrieval_sort(setdiff(1:numel(IDsort), nanID)));
        
        [firstadj2retrieval_sort, IDsort] = sort(firstadj2retrieval);
        IDsort(isnan(firstadj2retrieval_sort)) = [];
        subplot(2, 2, 3);
        nanID = heatmap4gcamp(t_aligned, zscore_retrievalstart_all(:, IDsort, i), 'Sort by hand adjustment');
        event4heatmap(t_aligned, firstadj2retrieval_sort(setdiff(1:numel(IDsort), nanID)));
        
        [firstbite2retrieval_sort, IDsort] = sort(firstbite2retrieval);
        IDsort(isnan(firstbite2retrieval_sort)) = [];
        subplot(2, 2, 4);
        nanID = heatmap4gcamp(t_aligned, zscore_retrievalstart_all(:, IDsort, i), 'Sort by bite');
        event4heatmap(t_aligned, firstbite2retrieval_sort(setdiff(1:numel(IDsort), nanID)));
    end
    
    if ~isempty(zscore_withdraw_all)
        figure('Name', ['aligned to Withdraw: ' lgdtext{i}]);
        subplot(2, 2, 1);
        heatmap4gcamp(t_aligned, zscore_withdraw_all(:, :, i), 'Sort by trial number');
        set(gca, 'ButtonDownFcn', @extract_figure);
        
        [firstadj2withdraw_sort, IDsort] = sort(firstadj2withdraw);
        IDsort(isnan(firstadj2withdraw_sort)) = [];
        subplot(2, 2, 2);
        nanID = heatmap4gcamp(t_aligned, zscore_withdraw_all(:, IDsort, i), 'Sort by hand adjustment');
        event4heatmap(t_aligned, firstadj2withdraw_sort(setdiff(1:numel(IDsort), nanID)));
        
        [firstbite2withdraw_sort, IDsort] = sort(firstbite2withdraw);
        IDsort(isnan(firstbite2withdraw_sort)) = [];
        subplot(2, 2, 3);
        nanID = heatmap4gcamp(t_aligned, zscore_withdraw_all(:, IDsort, i), 'Sort by bite');
        event4heatmap(t_aligned, firstbite2withdraw_sort(setdiff(1:numel(IDsort), nanID)));
        
        [firstchew2withdraw_sort, IDsort] = sort(firstchew2withdraw);
        IDsort(isnan(firstchew2withdraw_sort)) = [];
        subplot(2, 2, 4);
        nanID = heatmap4gcamp(t_aligned, zscore_withdraw_all(:, IDsort, i), 'Sort by chew');
        event4heatmap(t_aligned, firstchew2withdraw_sort(setdiff(1:numel(IDsort), nanID)));
    end
    
    if ~isempty(zscore_handadj_all)
        figure('Name', ['aligned to Hand Adjustment: ' lgdtext{i}]);
        subplot(1, 2, 1);
        heatmap4gcamp(t_aligned, zscore_handadj_all(:, :, i), 'Sort by trial number');
        set(gca, 'ButtonDownFcn', @extract_figure);
        
        [firstbite2adj_sort, IDsort] = sort(firstbite2adj);
        IDsort(isnan(firstbite2adj_sort)) = [];
        subplot(1, 2, 2);
        nanID = heatmap4gcamp(t_aligned, zscore_handadj_all(:, IDsort, i), 'Sort by bite');
        event4heatmap(t_aligned, firstbite2adj_sort(setdiff(1:numel(IDsort), nanID)));
    end
    
    if ~isempty(zscore_bite_all)
        figure('Name', ['aligned to Bite: ' lgdtext{i}]);
        subplot(1, 2, 1);
        heatmap4gcamp(t_aligned, zscore_bite_all(:, :, i), 'Sort by trial number');
        set(gca, 'ButtonDownFcn', @extract_figure);
        
        [firstadj2bite_sort, IDsort] = sort(firstadj2bite);
        IDsort(isnan(firstadj2bite_sort)) = [];
        subplot(1, 2, 2);
        nanID = heatmap4gcamp(t_aligned, zscore_bite_all(:, IDsort, i), 'Sort by hand adjustment');
        event4heatmap(t_aligned, firstadj2bite_sort(setdiff(1:numel(IDsort), nanID)));
    end
end

colors = [0.3010 0.7450 0.9330; 0.8500 0.3250 0.0980];
hp = zeros(1, nchannel);
if ~isempty(zscore_retrievalstart_all)
    figure;
    for i = 1:nchannel
        hp(i) = plot_tj_MeanSEM(t_aligned, zscore_retrievalstart_all(:, :, i), colors(i, :), colors(i, :), 'Time (s)', 'Z-score', 'aligned to Retrieval Start with Mouth');
    end
    yl = ylim;
    plot([0 0], yl, '--k', 'LineWidth', 1);
    plot(ones(1, 2)*mean(firstwithdraw2retrieval(firstwithdraw2retrieval <= post)), yl, '--k', 'LineWidth', 1);
    legend(hp, lgdtext, 'Location', 'NorthEast');
    legend('boxoff');
end

figure;
for i = 1:nchannel
%     if ~isempty(zscore_retrievalstart_all)
%         plot_tj_MeanSEM(t_aligned, zscore_retrievalstart_all(:, :, i), colors(i, :), colors(i, :), 'Time (s)', 'Z-score', 'aligned to All');
%     end
    if ~isempty(zscore_withdraw_all)
        plot_tj_MeanSEM(t_aligned+0*(pre+post)*(1+gap), zscore_withdraw_all(:, :, i), colors(i, :), colors(i, :), 'Time (s)', 'Z-score', 'aligned to All');
    end
    if ~isempty(zscore_handadj_all)
        plot_tj_MeanSEM(t_aligned+1*(pre+post)*(1+gap), zscore_handadj_all(:, :, i), colors(i, :), colors(i, :), 'Time (s)', 'Z-score', 'aligned to All');
    end
    if ~isempty(zscore_bite_all)
        plot_tj_MeanSEM(t_aligned+2*(pre+post)*(1+gap), zscore_bite_all(:, :, i), colors(i, :), colors(i, :), 'Time (s)', 'Z-score', 'aligned to All');
    end
    if ~isempty(zscore_chew_all)
        hp(i) = plot_tj_MeanSEM(t_aligned+3*(pre+post)*(1+gap), zscore_chew_all(:, :, i), colors(i, :), colors(i, :), 'Time (s)', 'Z-score', 'aligned to All');
    end
end
n = 4; % number of plots
plot([-pre post+(n-1)*(pre+post)*(1+gap)], [0 0], '--k', 'LineWidth', 1); % Z-score 0
xlim([-pre post+(n-1)*(pre+post)*(1+gap)]);
yl = ylim;
for i = 1:n
    plot([(i-1)*(pre+post)*(1+gap) (i-1)*(pre+post)*(1+gap)], yl, '--k', 'LineWidth', 1); % time 0
end
xtick = [-pre:(pre+post)*(1+gap):-pre+(n-1)*(pre+post)*(1+gap) 0:(pre+post)*(1+gap):(n-1)*(pre+post)*(1+gap) post:(pre+post)*(1+gap):post+(n-1)*(pre+post)*(1+gap)];
[xtick, ID] = sort(xtick);
xtick_text = repmat([-pre 0 post], n, 1);
xtick_text = xtick_text(ID);
set(gca, 'XTick', xtick, 'XTickLabel', num2str(xtick_text(:)));
legend(hp, lgdtext, 'Location', 'NorthEast');
legend('boxoff');

result.zscore_stepcross = zscore_stepcross_all;
result.zscore_retrievalstart = zscore_retrievalstart_all;
result.zscore_withdraw = zscore_withdraw_all;
result.zscore_handadj = zscore_handadj_all;
result.zscore_bite = zscore_bite_all;
result.zscore_chew = zscore_chew_all;
result.t = t_aligned;
result.retrievalstart2stepcross = retrievalstart2stepcross;
result.firstwithdraw2retrieval = firstwithdraw2retrieval;
result.firstadj2retrieval = firstadj2retrieval;
result.firstbite2retrieval = firstbite2retrieval;
result.firstadj2withdraw = firstadj2withdraw;
result.firstbite2withdraw = firstbite2withdraw;
result.firstchew2withdraw = firstchew2withdraw;
result.firstbite2adj = firstbite2adj;
result.firstadj2bite = firstadj2bite;
assignin('base', 'result', result);