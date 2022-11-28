function EventProperty(app, Exp_Path)
value = app.TrialsListBox.Value;
value = sort(value);
trials = numel(value);

adjustment2firstbite = [];
sit2retrievalend = [];
lefthandreach2retrievalend = [];
righthandreach2retrievalend = [];

firstadjustment2withdraw = [];
firstbite2withdraw = [];

sit2retrievalstart = [];
lefthandreach2retrievalstart = [];
righthandreach2retrievalstart = [];

feedingduration = [];

interbout_interval = [];

nadjustment_NI = 0;
nadjustment_I = 0;

% get bite events
Bite_events = [];
try
    audiolocation = Exp_Path(1:end-7);
    temp = load([audiolocation '\Detected_Bite_Events.mat']);
    Bite_events = temp.Audio_analysis;
end

try
    temp = load([Exp_Path '\Analysis_Session.mat']);
    video_annotation = temp.Video_annotation;
end

for i = 1:trials
    LabelledEvents = [];
    temp = load([Exp_Path '\LabelledEvents' num2str(value(i)) '.mat']);
    LabelledEvents = temp.LabelledEvents;
    
    if isfield(Bite_events(value(i)), 'laser_timestamps')
        if isempty(Bite_events(value(i)).laser_timestamps)
            nadjustment_NI = nadjustment_NI+sum(LabelledEvents.PawLAdjustmentStart >= 4 & LabelledEvents.PawLAdjustmentStart < 8)...
                +sum(LabelledEvents.PawRAdjustmentStart >= 4 & LabelledEvents.PawRAdjustmentStart < 8);
        else
            nadjustment_I = nadjustment_I+sum(LabelledEvents.PawLAdjustmentStart >= 4 & LabelledEvents.PawLAdjustmentStart < 8)...
                +sum(LabelledEvents.PawRAdjustmentStart >= 4 & LabelledEvents.PawRAdjustmentStart < 8);
        end
    end
    
    bite_timestamps = [];
    if ~isempty(Bite_events)
        bite_timestamps = Bite_events(value(i)).time_bites;
    end
    
    BiteBoutStart = [];
    if isfield(LabelledEvents, 'BiteBoutStartHMM') && ~isempty(LabelledEvents.BiteBoutStartHMM)
        BiteBoutStart = LabelledEvents.BiteBoutStartHMM;
    elseif ~isempty(LabelledEvents.BiteBoutStart)
        BiteBoutStart = LabelledEvents.BiteBoutStart;
    end
    if ~isempty(BiteBoutStart)
        BiteBoutStart = sort(BiteBoutStart);
        if numel(BiteBoutStart) > 1
            interbout_interval = [interbout_interval; diff(BiteBoutStart)];
        end
        PawRAdjustmentStart = LabelledEvents.PawRAdjustmentStart;
        PawLAdjustmentStart = LabelledEvents.PawLAdjustmentStart;
        for j = 1:numel(BiteBoutStart)
            currentbout = BiteBoutStart(j);
            if j ~= numel(BiteBoutStart)
                nextbout = BiteBoutStart(j+1);
            else
                nextbout = inf;
            end
            biteinbout = bite_timestamps(bite_timestamps > currentbout & bite_timestamps < nextbout);
            PawRAdjustmentinbout = PawRAdjustmentStart(PawRAdjustmentStart > currentbout & PawRAdjustmentStart < nextbout);
            PawLAdjustmentinbout = PawLAdjustmentStart(PawLAdjustmentStart > currentbout & PawLAdjustmentStart < nextbout);
            if ~isempty(biteinbout)
                firstbite2withdraw = [firstbite2withdraw; min(biteinbout)-currentbout];
                if ~isempty(PawRAdjustmentinbout)
                    adjustment2firstbite = [adjustment2firstbite; PawRAdjustmentinbout-min(biteinbout)];
                end
                if ~isempty(PawLAdjustmentinbout)
                    adjustment2firstbite = [adjustment2firstbite; PawLAdjustmentinbout-min(biteinbout)];
                end
            end
            
            if ~(isempty(PawRAdjustmentinbout) && isempty(PawLAdjustmentinbout))
                adjustment2withdraw = [min(PawRAdjustmentinbout)-currentbout min(PawLAdjustmentinbout)-currentbout];
                firstadjustment2withdraw = [firstadjustment2withdraw; min(adjustment2withdraw)];
            end
        end
    end
    
    RetrievalStart = LabelledEvents.RetrievalStart;
    if ~isempty(RetrievalStart)
        RetrievalStart = sort(RetrievalStart);
        MouthRetrievalStart = LabelledEvents.MouthRetrievalStart;
        MouthRetrievalEnd = LabelledEvents.MouthRetrievalEnd;
        SitStart = LabelledEvents.SitStart;
        SitStart = sort(SitStart);
        PawLReachStart = LabelledEvents.PawLReachStart;
        PawLReachStart = sort(PawLReachStart);
        PawRReachStart = LabelledEvents.PawRReachStart;
        PawRReachStart = sort(PawRReachStart);
        if (numel(SitStart) == numel(PawLReachStart)) && (numel(SitStart) == numel(PawRReachStart))
            if numel(SitStart) < numel(RetrievalStart)
                IDmatch = nan(1, numel(SitStart));
                for j = 1:numel(SitStart)
                    [~, IDmatch(j)] = min(abs(RetrievalStart-SitStart(j)));
                end
                RetrievalStart = RetrievalStart(IDmatch);
            end
            
            for j = 1:numel(RetrievalStart)
                MouthRetrievalID = find(MouthRetrievalStart == RetrievalStart(j));
                if ~isempty(MouthRetrievalID)
                    [~, SitID] = min(abs(SitStart-MouthRetrievalEnd(MouthRetrievalID)));
                    sit2retrievalend = [sit2retrievalend; SitStart(SitID)-MouthRetrievalEnd(MouthRetrievalID)];
                    lefthandreach2retrievalend = [lefthandreach2retrievalend; PawLReachStart(SitID)-MouthRetrievalEnd(MouthRetrievalID)];
                    righthandreach2retrievalend = [righthandreach2retrievalend; PawRReachStart(SitID)-MouthRetrievalEnd(MouthRetrievalID)];
                    
                    sit2retrievalstart = [sit2retrievalstart; SitStart(SitID)-RetrievalStart(j)];
                    lefthandreach2retrievalstart = [lefthandreach2retrievalstart; PawLReachStart(SitID)-RetrievalStart(j)];
                    righthandreach2retrievalstart = [righthandreach2retrievalstart; PawRReachStart(SitID)-RetrievalStart(j)];
                end
            end
        end
    end
    
    FeedingEnd = LabelledEvents.FeedingEnd;
    if ~isempty(FeedingEnd)
        feedingduration = [feedingduration; FeedingEnd];
    elseif isempty(FeedingEnd)
        feedingduration = [feedingduration; video_annotation(value(i)).time_feeding_end];
    end
end
result.adjustment2firstbite = adjustment2firstbite;
result.sit2retrievalend = sit2retrievalend;
result.lefthandreach2retrievalend = lefthandreach2retrievalend;
result.righthandreach2retrievalend = righthandreach2retrievalend;
result.firstbite2withdraw = firstbite2withdraw;
result.firstadjustment2withdraw = firstadjustment2withdraw;
result.sit2retrievalstart = sit2retrievalstart;
result.lefthandreach2retrievalstart = lefthandreach2retrievalstart;
result.righthandreach2retrievalstart = righthandreach2retrievalstart;
result.feedingduration = feedingduration;
result.interbout_interval = interbout_interval;

result.nadjustment_NI = nadjustment_NI;
result.nadjustment_I = nadjustment_I;
assignin('base', 'result', result);
msgbox('Done !');