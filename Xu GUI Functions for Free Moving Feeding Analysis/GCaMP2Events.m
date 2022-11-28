function GCaMP2Events(app, Exp_Path)
value = app.TrialsListBox.Value;
value = sort(value);
trials = numel(value);

pre = app.presEditField.Value;
post = app.postsEditField.Value;

zscore_tongueout_all = [];
zscore_tongueout_random_all = [];
zscore_pawlreach_all = [];
zscore_pawlreach_random_all = [];
zscore_pawrreach_all = [];
zscore_pawrreach_random_all = [];
zscore_sit_all = [];
zscore_sit_random_all = [];
zscore_pawladjustment_all = [];
zscore_pawladjustment_random_all = [];
zscore_pawradjustment_all = [];
zscore_pawradjustment_random_all = [];

zscore_retrievalstart_all = [];
zscore_retrievalstart_random_all = [];
zscore_biteboutstart_all = [];
zscore_biteboutstart_random_all = [];
zscore_biteboutstartAdj_all = [];
zscore_biteboutstartAdj_random_all = [];
zscore_biteboutstartNAdj_all = [];
zscore_biteboutstartNAdj_random_all = [];
zscore_lastbite_all = []; % chew start
zscore_lastbite_random_all = [];
zscore_firstbite_all = [];
zscore_firstbite_random_all = [];
zscore_barriercross_all = []; % time 0
zscore_barriercross_random_all = [];

RetrievalStartMouth_all = [];
SitEnd2RetrievalStartMouth_all = [];
lastbite2firstbite_all = [];
firstbite2lastbite_all = [];
biteboutstart2firstbite_all = [];
biteboutstart2lastbite_all = [];
lastbite2biteboutstart_all = [];

zscore_biteboutstart_groups = {[] [] []};
ngroup = numel(zscore_biteboutstart_groups);

newtimestamps = [];

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
%             %
%             if isfield(events, 'BiteBoutStartHMM')
%                 events.BiteBoutStartHMM = [];
%             end
%             %
        catch
            continue;
        end
        
        bite_timestamps = [];
        if ~isempty(Bite_events)
            bite_timestamps = Bite_events(value(i)).time_bites;
        end
        
        fpdata = fpdata_all.zsignal_all(:, value(i));
        fpdata_zsignal = [];
        fpdata_t = [];
        for j = 1:nchannel
            fpdata_zsignal(:, j) = fpdata{j}(:, 2);
        end
        fpdata_t = fpdata{1}(:, 1);
        SampleRate = mean(diff(fpdata_t));
        
        % these are not in use
        if ~isempty(events.TongueOut)
            [t_aligned, zscore_tongueout, zscore_tongueout_random] = AlignSignal2Event(fpdata_t, fpdata_zsignal, events.TongueOut, pre, post);
            zscore_tongueout_all = [zscore_tongueout_all zscore_tongueout];
            zscore_tongueout_random_all = [zscore_tongueout_random_all zscore_tongueout_random];
        end
        if ~isempty(events.PawLReachStart)
            [t_aligned, zscore_pawlreach, zscore_pawlreach_random] = AlignSignal2Event(fpdata_t, fpdata_zsignal, events.PawLReachStart, pre, post);
            zscore_pawlreach_all = [zscore_pawlreach_all zscore_pawlreach];
            zscore_pawlreach_random_all = [zscore_pawlreach_random_all zscore_pawlreach_random];
        end
        if ~isempty(events.PawRReachStart)
            [t_aligned, zscore_pawrreach, zscore_pawrreach_random] = AlignSignal2Event(fpdata_t, fpdata_zsignal, events.PawRReachStart, pre, post);
            zscore_pawrreach_all = [zscore_pawrreach_all zscore_pawrreach];
            zscore_pawrreach_random_all = [zscore_pawrreach_random_all zscore_pawrreach_random];
        end
        if ~isempty(events.SitStart)
            [t_aligned, zscore_sit, zscore_sit_random] = AlignSignal2Event(fpdata_t, fpdata_zsignal, events.SitStart, pre, post);
            zscore_sit_all = [zscore_sit_all zscore_sit];
            zscore_sit_random_all = [zscore_sit_random_all zscore_sit_random];
        end
        
        % possible useful ones start from here
        if ~isempty(events.PawLAdjustmentStart)
            newtimestamps = get_ramdomized_time_in_handling_stage(events, bite_timestamps);
            if ~isempty(newtimestamps.PawLAdjustmentStart)
                [t_aligned, zscore_pawladjustment, ~] = AlignSignal2Event(fpdata_t, fpdata_zsignal, newtimestamps.PawLAdjustmentStart, pre, post);
                [t_aligned, zscore_pawladjustment_random, ~] = AlignSignal2Event(fpdata_t, fpdata_zsignal, newtimestamps.PawLAdjustmentRnd, pre, post);
                zscore_pawladjustment_all = [zscore_pawladjustment_all zscore_pawladjustment];
                zscore_pawladjustment_random_all = [zscore_pawladjustment_random_all zscore_pawladjustment_random];
            end
        end
        if ~isempty(events.PawRAdjustmentStart)
            if isempty(newtimestamps)
                newtimestamps = get_ramdomized_time_in_handling_stage(events, bite_timestamps);
            end
            if ~isempty(newtimestamps.PawRAdjustmentStart)
                [t_aligned, zscore_pawradjustment, ~] = AlignSignal2Event(fpdata_t, fpdata_zsignal, newtimestamps.PawRAdjustmentStart, pre, post);
                [t_aligned, zscore_pawradjustment_random, ~] = AlignSignal2Event(fpdata_t, fpdata_zsignal, newtimestamps.PawRAdjustmentRnd, pre, post);
                zscore_pawradjustment_all = [zscore_pawradjustment_all zscore_pawradjustment];
                zscore_pawradjustment_random_all = [zscore_pawradjustment_random_all zscore_pawradjustment_random];
            end
        end
        
        if ~isempty(events.RetrievalStart)
            if numel(events.SitEnd) == numel(events.RetrievalStart)
                events.SitEnd = events.SitEnd(ismember(events.RetrievalStart, events.MouthRetrievalStart));
            else
                events.SitEnd = [];
            end
            events.RetrievalStart = events.RetrievalStart(ismember(events.RetrievalStart, events.MouthRetrievalStart)); % only include successful retrieval using mouth
            if ~isempty(events.SitEnd)
                SitEnd2RetrievalStartMouth_all = [SitEnd2RetrievalStartMouth_all; events.SitEnd-events.RetrievalStart];
            end
            RetrievalStartMouth_all = [RetrievalStartMouth_all; events.RetrievalStart];
            if ~isempty(events.RetrievalStart)
                [t_aligned, zscore_retrievalstart, zscore_retrievalstart_random] = AlignSignal2Event(fpdata_t, fpdata_zsignal, events.RetrievalStart, pre, post);
                zscore_retrievalstart_all = [zscore_retrievalstart_all zscore_retrievalstart];
                zscore_retrievalstart_random_all = [zscore_retrievalstart_random_all zscore_retrievalstart_random];
            end
        end
        
        if isfield(events, 'BiteBoutStartHMM') && ~isempty(events.BiteBoutStartHMM)
            [t_aligned, zscore_biteboutstart, zscore_biteboutstart_random] = AlignSignal2Event(fpdata_t, fpdata_zsignal, events.BiteBoutStartHMM, pre, post);
            zscore_biteboutstart_all = [zscore_biteboutstart_all zscore_biteboutstart];
            zscore_biteboutstart_random_all = [zscore_biteboutstart_random_all zscore_biteboutstart_random];
        elseif ~isempty(events.BiteBoutStart)
            [t_aligned, zscore_biteboutstart, zscore_biteboutstart_random] = AlignSignal2Event(fpdata_t, fpdata_zsignal, events.BiteBoutStart, pre, post);
            zscore_biteboutstart_all = [zscore_biteboutstart_all zscore_biteboutstart];
            zscore_biteboutstart_random_all = [zscore_biteboutstart_random_all zscore_biteboutstart_random];
        end
        
        % separate the times in different groups
        if isfield(events, 'BiteBoutStartHMM') && ~isempty(events.BiteBoutStartHMM)
            nBiteBoutStart = numel(events.BiteBoutStartHMM);
            events.BiteBoutStartHMM = sort(events.BiteBoutStartHMM);
            if nBiteBoutStart >= ngroup
                for j = 1:ngroup
                    [t_aligned, zscore_biteboutstart, ~] = AlignSignal2Event(fpdata_t, fpdata_zsignal,...
                        events.BiteBoutStartHMM((j-1)*floor(nBiteBoutStart/ngroup)+1:j*floor(nBiteBoutStart/ngroup)), pre, post);
                    zscore_biteboutstart_groups{j} = [zscore_biteboutstart_groups{j} zscore_biteboutstart];
                end
            end
        elseif ~isempty(events.BiteBoutStart)
            nBiteBoutStart = numel(events.BiteBoutStart);
            events.BiteBoutStart = sort(events.BiteBoutStart);
            if nBiteBoutStart >= ngroup
                for j = 1:ngroup
                    [t_aligned, zscore_biteboutstart, ~] = AlignSignal2Event(fpdata_t, fpdata_zsignal,...
                        events.BiteBoutStart((j-1)*floor(nBiteBoutStart/ngroup)+1:j*floor(nBiteBoutStart/ngroup)), pre, post);
                    zscore_biteboutstart_groups{j} = [zscore_biteboutstart_groups{j} zscore_biteboutstart];
                end
            end
        end
        
        % group the times based on whether there's paw adjustment in the bout
        if ~isempty(events.PawLAdjustmentStart)
            if ~isempty(events.FeedingEnd)
                bitebout = [events.RetrievalStart; events.BiteBoutStart];
            else
                bitebout = [events.RetrievalStart; events.BiteBoutStart(1:end-1)];
            end
            bitebout = sort(bitebout);
            RetrievalStartID = find(ismember(bitebout, events.RetrievalStart));
            AdjBoutID = true(size(bitebout));
            for j = 1:numel(bitebout)
                if j ~= numel(bitebout)
                    checkAdj = events.PawLAdjustmentStart > bitebout(j) & events.PawLAdjustmentStart < bitebout(j+1);
                else
                    checkAdj = events.PawLAdjustmentStart > bitebout(j);
                end
                if ~any(checkAdj)
                    AdjBoutID(j) = false;
                end
            end
            AdjBoutID(RetrievalStartID) = [];
            if ~isempty(events.FeedingEnd)
                if any(AdjBoutID)
                    [t_aligned, zscore_biteboutstart, zscore_biteboutstart_random] = AlignSignal2Event(fpdata_t, fpdata_zsignal, events.BiteBoutStart(AdjBoutID), pre, post);
                    zscore_biteboutstartAdj_all = [zscore_biteboutstartAdj_all zscore_biteboutstart];
                    zscore_biteboutstartAdj_random_all = [zscore_biteboutstartAdj_random_all zscore_biteboutstart_random];
                end
                if any(~AdjBoutID)
                    [t_aligned, zscore_biteboutstart, zscore_biteboutstart_random] = AlignSignal2Event(fpdata_t, fpdata_zsignal, events.BiteBoutStart(~AdjBoutID), pre, post);
                    zscore_biteboutstartNAdj_all = [zscore_biteboutstartNAdj_all zscore_biteboutstart];
                    zscore_biteboutstartNAdj_random_all = [zscore_biteboutstartNAdj_random_all zscore_biteboutstart_random];
                end
            else
                if any(AdjBoutID)
                    [t_aligned, zscore_biteboutstart, zscore_biteboutstart_random] = AlignSignal2Event(fpdata_t, fpdata_zsignal, events.BiteBoutStart([AdjBoutID; false]), pre, post);
                    zscore_biteboutstartAdj_all = [zscore_biteboutstartAdj_all zscore_biteboutstart];
                    zscore_biteboutstartAdj_random_all = [zscore_biteboutstartAdj_random_all zscore_biteboutstart_random];
                end
                if any(~AdjBoutID)
                    [t_aligned, zscore_biteboutstart, zscore_biteboutstart_random] = AlignSignal2Event(fpdata_t, fpdata_zsignal, events.BiteBoutStart([~AdjBoutID; false]), pre, post);
                    zscore_biteboutstartNAdj_all = [zscore_biteboutstartNAdj_all zscore_biteboutstart];
                    zscore_biteboutstartNAdj_random_all = [zscore_biteboutstartNAdj_random_all zscore_biteboutstart_random];
                end
            end
        end
        
        bitebout = [];
        if ~isempty(bite_timestamps)
            if ~isempty(events.FeedingEnd) || (isfield(events, 'BiteBoutStartHMM') && ~isempty(events.BiteBoutStartHMM))
                if isfield(events, 'BiteBoutStartHMM') && ~isempty(events.BiteBoutStartHMM)
                    bitebout = [events.RetrievalStart; events.BiteBoutStartHMM; inf];
                else
                    bitebout = [events.RetrievalStart; events.BiteBoutStart; inf];
                end
            elseif ~isempty(events.BiteBoutStart)
                bitebout = [events.RetrievalStart; events.BiteBoutStart];
            end
            if ~isempty(bitebout)
                bitebout = sort(bitebout);
                RetrievalStartID = find(ismember(bitebout, events.RetrievalStart));
                ninterval = numel(bitebout)-1;
                lastbite_timestamps = nan(ninterval, 1);
                firstbite_timestamps = nan(ninterval, 1);
                for k = 1:ninterval
                    lastbiteID = find(bite_timestamps > bitebout(k) & bite_timestamps < bitebout(k+1), 1, 'last');
                    firstbiteID = find(bite_timestamps > bitebout(k) & bite_timestamps < bitebout(k+1), 1, 'first');
                    if ~isempty(lastbiteID)
                        lastbite_timestamps(k) = bite_timestamps(lastbiteID);
                        firstbite_timestamps(k) = bite_timestamps(firstbiteID);
                        
                        if ~ismember(k, RetrievalStartID)
                            biteboutstart2firstbite_all = [biteboutstart2firstbite_all firstbite_timestamps(k)-bitebout(k)];
                            biteboutstart2lastbite_all = [biteboutstart2lastbite_all lastbite_timestamps(k)-bitebout(k)];
                        end
                        if k ~= ninterval && ~ismember(k+1, RetrievalStartID)
                            lastbite2biteboutstart_all = [lastbite2biteboutstart_all bitebout(k+1)-lastbite_timestamps(k)];
                        end
                    end
                end
                firstbite2lastbite_all = lastbite_timestamps-firstbite_timestamps;
                if max(RetrievalStartID) > ninterval
                    firstbite2lastbite_all(RetrievalStartID(1:end-1)) = []; % ignore the very first bite after each retrieval
                else
                    firstbite2lastbite_all(RetrievalStartID) = []; % ignore the very first bite after each retrieval
                end
                firstbite2lastbite_all(isnan(firstbite2lastbite_all)) = [];
                lastbite2firstbite_all = firstbite_timestamps(2:end)-lastbite_timestamps(1:end-1);
                if numel(RetrievalStartID) > 1
                    if max(RetrievalStartID) > ninterval
                        lastbite2firstbite_all(RetrievalStartID(2:end-1)-1) = [];
                    else
                        lastbite2firstbite_all(RetrievalStartID(2:end)-1) = [];
                    end
                end
                lastbite2firstbite_all(isnan(lastbite2firstbite_all)) = [];
                if numel(RetrievalStartID) > 1
                    if max(RetrievalStartID) > ninterval
                        lastbite_timestamps(RetrievalStartID(2:end-1)-1) = [];
                    else
                        lastbite_timestamps(RetrievalStartID(2:end)-1) = [];
                    end
                end
                lastbite_timestamps(isnan(lastbite_timestamps)) = [];
                if max(RetrievalStartID) > ninterval
                    firstbite_timestamps(RetrievalStartID(1:end-1)) = []; % discard the very first bite after each retrieval
                else
                    firstbite_timestamps(RetrievalStartID) = []; % discard the very first bite after each retrieval
                end
                firstbite_timestamps(isnan(firstbite_timestamps)) = [];
                if ~isempty(lastbite_timestamps)
                    [t_aligned, zscore_lastbite, zscore_lastbite_random] = AlignSignal2Event(fpdata_t, fpdata_zsignal, lastbite_timestamps, pre, post);
                    zscore_lastbite_all = [zscore_lastbite_all zscore_lastbite];
                    zscore_lastbite_random_all = [zscore_lastbite_random_all zscore_lastbite_random];
                end
                if ~isempty(firstbite_timestamps)
                    [t_aligned, zscore_firstbite, zscore_firstbite_random] = AlignSignal2Event(fpdata_t, fpdata_zsignal, firstbite_timestamps, pre, post);
                    zscore_firstbite_all = [zscore_firstbite_all zscore_firstbite];
                    zscore_firstbite_random_all = [zscore_firstbite_random_all zscore_firstbite_random];
                end
            end
        end

        % aligned to time 0 for barrier cross
        [t_aligned, zscore_barriercross, zscore_barriercross_random] = AlignSignal2Event(fpdata_t, fpdata_zsignal, 0, pre, post);
        zscore_barriercross_all = [zscore_barriercross_all zscore_barriercross];
        zscore_barriercross_random_all = [zscore_barriercross_random_all zscore_barriercross_random];
    end
end

for i = 1:nchannel
    figure;
    if ~isempty(zscore_tongueout_all)
        subplot(3, 3, 1);
        plot_gcamp_align2event(t_aligned, zscore_tongueout_random_all(:, :, i), zscore_tongueout_all(:, :, i), ['aligned to Tongue Out: ' lgdtext{i}]);
        set(gca, 'ButtonDownFcn', @extract_figure);
    end
    
    if ~isempty(zscore_pawlreach_all)
        subplot(3, 3, 2);
        plot_gcamp_align2event(t_aligned, zscore_pawlreach_random_all(:, :, i), zscore_pawlreach_all(:, :, i), ['aligned to Left Paw Reach: ' lgdtext{i}]);
        set(gca, 'ButtonDownFcn', @extract_figure);
    end
    
    if ~isempty(zscore_pawrreach_all)
        subplot(3, 3, 3);
        plot_gcamp_align2event(t_aligned, zscore_pawrreach_random_all(:, :, i), zscore_pawrreach_all(:, :, i), ['aligned to Right Paw Reach: ' lgdtext{i}]);
        set(gca, 'ButtonDownFcn', @extract_figure);
    end
    
    if ~isempty(zscore_sit_all)
        subplot(3, 3, 4);
        plot_gcamp_align2event(t_aligned, zscore_sit_random_all(:, :, i), zscore_sit_all(:, :, i), ['aligned to Sit: ' lgdtext{i}]);
        set(gca, 'ButtonDownFcn', @extract_figure);
    end
    
    if ~isempty(zscore_pawladjustment_all)
        subplot(3, 3, 5);
        plot_gcamp_align2event(t_aligned, zscore_pawladjustment_random_all(:, :, i), zscore_pawladjustment_all(:, :, i), ['aligned to Left Paw Adjustment: ' lgdtext{i}]);
        set(gca, 'ButtonDownFcn', @extract_figure);
    end
    
    if ~isempty(zscore_pawradjustment_all)
        subplot(3, 3, 6);
        plot_gcamp_align2event(t_aligned, zscore_pawradjustment_random_all(:, :, i), zscore_pawradjustment_all(:, :, i), ['aligned to Right Paw Adjustment: ' lgdtext{i}]);
        set(gca, 'ButtonDownFcn', @extract_figure);
    end
    
    if ~isempty(zscore_retrievalstart_all)
        subplot(3, 3, 7);
        plot_gcamp_align2event(t_aligned, zscore_retrievalstart_random_all(:, :, i), zscore_retrievalstart_all(:, :, i), ['aligned to Retrieval Start with Mouth: ' lgdtext{i}]);
        set(gca, 'ButtonDownFcn', @extract_figure);
    end
    
    if ~isempty(zscore_biteboutstart_all)
        subplot(3, 3, 8);
        plot_gcamp_align2event(t_aligned, zscore_biteboutstart_random_all(:, :, i), zscore_biteboutstart_all(:, :, i), ['aligned to Bite Bout Start: ' lgdtext{i}]);
        set(gca, 'ButtonDownFcn', @extract_figure);
    end
    
    if ~isempty(zscore_lastbite_all)
        subplot(3, 3, 9);
        plot_gcamp_align2event(t_aligned, zscore_lastbite_random_all(:, :, i), zscore_lastbite_all(:, :, i), ['aligned to last bite of each bite bout: ' lgdtext{i}]);
        set(gca, 'ButtonDownFcn', @extract_figure);
    end
end

for i = 1:nchannel
    if ~isempty(zscore_biteboutstart_all)
        figure;
        heatmap4gcamp(t_aligned, zscore_biteboutstart_all(:, :, i), ['aligned to Bite Bout Start: ' lgdtext{i}]);
    end
    
    if ~isempty(zscore_lastbite_all)
        figure;
        heatmap4gcamp(t_aligned, zscore_lastbite_all(:, :, i), ['aligned to last bite of each bite bout: ' lgdtext{i}]);
        
        figure;
        heatmap4gcamp(t_aligned, zscore_firstbite_all(:, :, i), ['aligned to first bite of each bite bout: ' lgdtext{i}]);
    end
    
    if ~isempty(zscore_retrievalstart_all)
        figure;
        heatmap4gcamp(t_aligned, zscore_retrievalstart_all(:, :, i), ['aligned to Retrieval Start with Mouth: ' lgdtext{i}]);
    end
    
    figure;
    heatmap4gcamp(t_aligned, zscore_barriercross_all(:, :, i), ['aligned to Barrier Cross: ' lgdtext{i}]);
end

colors = [0.3010 0.7450 0.9330; 0.8500 0.3250 0.0980];
hp = zeros(1, nchannel);
figure;
if ~isempty(zscore_biteboutstart_all)
    subplot(3, 3, 1);
    for i = 1:nchannel
        hp(i) = plot_tj_MeanSEM(t_aligned, zscore_biteboutstart_all(:, :, i), colors(i, :), colors(i, :), 'Time (s)', 'Z-score', 'aligned to Bite Bout Start');
    end
    yl = ylim;
    plot([0 0], yl, '--k', 'LineWidth', 1);
    if ~isempty(biteboutstart2firstbite_all)
        plot(ones(1, 2)*mean(biteboutstart2firstbite_all(biteboutstart2firstbite_all <= post)), yl, '--k', 'LineWidth', 1);
        plot(ones(1, 2)*mean(biteboutstart2lastbite_all(biteboutstart2lastbite_all <= post)), yl, '--k', 'LineWidth', 1);
    end
    legend(hp, lgdtext, 'Location', 'NorthEast');
    legend('boxoff');
end

if ~isempty(zscore_lastbite_all)
    subplot(3, 3, 2);
    for i = 1:nchannel
        hp(i) = plot_tj_MeanSEM(t_aligned, zscore_lastbite_all(:, :, i), colors(i, :), colors(i, :), 'Time (s)', 'Z-score', 'aligned to last bite of each bite bout');
    end
    yl = ylim;
    plot([0 0], yl, '--k', 'LineWidth', 1);
    plot(ones(1, 2)*mean(lastbite2biteboutstart_all(lastbite2biteboutstart_all <= post)), yl, '--k', 'LineWidth', 1);
    plot(ones(1, 2)*mean(lastbite2firstbite_all(lastbite2firstbite_all <= post)), yl, '--k', 'LineWidth', 1);
    legend(hp, lgdtext, 'Location', 'NorthEast');
    legend('boxoff');
    
    subplot(3, 3, 3);
    for i = 1:nchannel
        hp(i) = plot_tj_MeanSEM(t_aligned, zscore_firstbite_all(:, :, i), colors(i, :), colors(i, :), 'Time (s)', 'Z-score', 'aligned to first bite of each bite bout');
    end
    yl = ylim;
    plot([0 0], yl, '--k', 'LineWidth', 1);
    plot(ones(1, 2)*mean(firstbite2lastbite_all(firstbite2lastbite_all <= post)), yl, '--k', 'LineWidth', 1);
    legend(hp, lgdtext, 'Location', 'NorthEast');
    legend('boxoff');
end

if ~isempty(zscore_retrievalstart_all)
    subplot(3, 3, 4);
    for i = 1:nchannel
        hp(i) = plot_tj_MeanSEM(t_aligned, zscore_retrievalstart_all(:, :, i), colors(i, :), colors(i, :), 'Time (s)', 'Z-score', 'aligned to Retrieval Start with Mouth');
    end
    yl = ylim;
    plot([0 0], yl, '--k', 'LineWidth', 1);
    plot(ones(1, 2)*mean(SitEnd2RetrievalStartMouth_all(SitEnd2RetrievalStartMouth_all <= post)), yl, '--k', 'LineWidth', 1);
    legend(hp, lgdtext, 'Location', 'NorthEast');
    legend('boxoff');
end

subplot(3, 3, 5);
for i = 1:nchannel
    hp(i) = plot_tj_MeanSEM(t_aligned, zscore_barriercross_all(:, :, i), colors(i, :), colors(i, :), 'Time (s)', 'Z-score', 'aligned to Barrier Cross');
end
yl = ylim;
plot([0 0], yl, '--k', 'LineWidth', 1);
plot(ones(1, 2)*mean(RetrievalStartMouth_all(RetrievalStartMouth_all <= post)), yl, '--k', 'LineWidth', 1);
legend(hp, lgdtext, 'Location', 'NorthEast');
legend('boxoff');

colors = [0 0 1; 0.3010 0.7450 0.9330; 0 1 1];
hp = zeros(1, ngroup);
for i = 1:nchannel
    figure('Name', lgdtext{i});
    for j = 1:ngroup
        hp(j) = plot_tj_MeanSEM(t_aligned, zscore_biteboutstart_groups{j}(:, :, i), colors(j, :), colors(j, :), 'Time (s)', 'Z-score', 'aligned to Bite Bout Start');
    end
    yl = ylim;
    plot([0 0], yl, '--k', 'LineWidth', 1);
    legend(hp, {'First', 'Middle', 'Last'}, 'Location', 'NorthEast');
    legend('boxoff');
end

if ~isempty(zscore_biteboutstartAdj_all) && ~isempty(zscore_biteboutstartNAdj_all)
    for i = 1:nchannel
        figure('Name', lgdtext{i});
        hp(1) = plot_tj_MeanSEM(t_aligned, zscore_biteboutstartAdj_all(:, :, i), colors(1, :), colors(1, :), 'Time (s)', 'Z-score', 'aligned to Bite Bout Start');
        hp(2) = plot_tj_MeanSEM(t_aligned, zscore_biteboutstartNAdj_all(:, :, i), colors(2, :), colors(2, :), 'Time (s)', 'Z-score', 'aligned to Bite Bout Start');
        yl = ylim;
        plot([0 0], yl, '--k', 'LineWidth', 1);
        legend(hp, {'with Adj', 'without Adj'}, 'Location', 'NorthEast');
        legend('boxoff');
    end
end

result.zscore_retrievalstart = zscore_retrievalstart_all;
result.zscore_biteboutstart = zscore_biteboutstart_all;
result.zscore_biteboutstartAdj = zscore_biteboutstartAdj_all;
result.zscore_biteboutstartNAdj = zscore_biteboutstartNAdj_all;
result.zscore_barriercross = zscore_barriercross_all;
result.zscore_lastbite = zscore_lastbite_all;
result.t = t_aligned;
assignin('base', 'result', result);
assignin('base', 'SitEnd2RetrievalStartMouth', SitEnd2RetrievalStartMouth_all);