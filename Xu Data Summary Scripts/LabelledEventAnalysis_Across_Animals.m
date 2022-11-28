%*---------------------------------------------------------------------*
% Huang Lab
% Cold Spring Harbor Laboratory
% Author : Xu An, Apr 2021
% xan@cshl.edu
% Version: 1.0
%
% Update 1
% add options for different paradigms by Xu An, June 2021
% Huang Lab
% Duke University
% xu.an@duke.edu
% Version: 2.0
%*---------------------------------------------------------------------*
%%
paradigm = 2; % 1 for 4-8s inhibition; 2 for 0-4s inhibition
[filename, pathname] = uigetfile('*.mat', 'Pick all of the individual data set', 'MultiSelect', 'on');
if isequal(filename, 0)
    return;
end
clc;
N = numel(filename);
APawLAdjustment_I = zeros(1, N); % 'A' means average
APawLAdjustment_NI = zeros(1, N);
APawRAdjustment_I = zeros(1, N);
APawRAdjustment_NI = zeros(1, N);
APawLRAdjustment_I = zeros(1, N);
APawLRAdjustment_NI = zeros(1, N);

AllPawLAdjustment_I = [];
AllPawLAdjustment_NI = [];
AllPawRAdjustment_I = [];
AllPawRAdjustment_NI = [];
AllPawLRAdjustment_I = [];
AllPawLRAdjustment_NI = [];
switch paradigm
    case 1
        ABiteBoutStart_I = zeros(1, N);
        ABiteBoutStart_NI = zeros(1, N);
        
        alladj_NI = zeros(1, N);
        alladj_I = zeros(1, N);
        biadj_NI = zeros(1, N);
        biadj_I = zeros(1, N);
        nbite_NI = zeros(1, N);
        nbite_I = zeros(1, N);
    case 2
        AMouthMovement_NI = zeros(1, N);
        AMouthMovement_I = zeros(1, N);
        ARetrievalStartT_NI = zeros(1, N);
        ARetrievalStartT_I = zeros(1, N);
        ARetrievalEndT_NI = zeros(1, N);
        ARetrievalEndT_I = zeros(1, N);
        AMouthRetrievalStart_NI = zeros(1, N);
        AMouthRetrievalStart_I = zeros(1, N);
        ASitDelayT_NI = zeros(1, N);
        ASitDelayT_I = zeros(1, N);
        APawLReachDelayT_NI = zeros(1, N);
        APawLReachDelayT_I = zeros(1, N);
        APawRReachDelayT_NI = zeros(1, N);
        APawRReachDelayT_I = zeros(1, N);
        ASitT_NI = zeros(1, N);
        ASitT_I = zeros(1, N);
        APawLReachT_NI = zeros(1, N);
        APawLReachT_I = zeros(1, N);
        APawRReachT_NI = zeros(1, N);
        APawRReachT_I = zeros(1, N);
        ABiteDelayT_NI = zeros(1, N);
        ABiteDelayT_I = zeros(1, N);
        ASitEndT1_NI = zeros(1, N);
        ASitEndT1_I = zeros(1, N);
end

for i = 1:N
    PawLAdjustment_NI = [];
    PawLAdjustment_I = [];
    PawRAdjustment_NI = [];
    PawRAdjustment_I = [];
    PawLRAdjustment_NI = [];
    PawLRAdjustment_I = [];
    switch paradigm
        case 1
            BiteBoutStart_I = [];
            BiteBoutStart_NI = [];
        case 2
            MouthMovement_NI = [];
            MouthMovement_I = [];
            RetrievalStartT_NI = [];
            RetrievalStartT_I = [];
            RetrievalEndT_NI = [];
            RetrievalEndT_I = [];
            MouthRetrievalStart_NI = [];
            MouthRetrievalStart_I = [];
            SitDelayT_NI = [];
            SitDelayT_I = [];
            SitT_NI = [];
            SitT_I = [];
            PawLReachDelayT_NI = [];
            PawLReachDelayT_I = [];
            PawLReachT_NI = [];
            PawLReachT_I = [];
            PawRReachDelayT_NI = [];
            PawRReachDelayT_I = [];
            PawRReachT_NI = [];
            PawRReachT_I = [];
            BiteDelayT_NI = [];
            BiteDelayT_I = [];
            SitEndT1_NI = [];
            SitEndT1_I = [];
    end
    load([pathname filename{i}]);
    for j = 1:numel(trials_NI)
        MouthOpen = trials_NI(j).LabelledEvents.MouthOpen;
        TongueOut = trials_NI(j).LabelledEvents.TongueOut;
        FoodinMouth = trials_NI(j).LabelledEvents.FoodinMouth;
        SitStart = trials_NI(j).LabelledEvents.SitStart;
        SitEnd = trials_NI(j).LabelledEvents.SitEnd;
        PawLReachStart = trials_NI(j).LabelledEvents.PawLReachStart;
        PawLReachEnd = trials_NI(j).LabelledEvents.PawLReachEnd;
        PawRReachStart = trials_NI(j).LabelledEvents.PawRReachStart;
        PawRReachEnd = trials_NI(j).LabelledEvents.PawRReachEnd;
        MouthRetrievalStart = trials_NI(j).LabelledEvents.MouthRetrievalStart;
        MouthRetrievalEnd = trials_NI(j).LabelledEvents.MouthRetrievalEnd;
        RetrievalStart = trials_NI(j).LabelledEvents.RetrievalStart;
        PawLAdjustmentStart = trials_NI(j).LabelledEvents.PawLAdjustmentStart;
        PawLAdjustmentEnd = trials_NI(j).LabelledEvents.PawLAdjustmentEnd;
        PawRAdjustmentStart = trials_NI(j).LabelledEvents.PawRAdjustmentStart;
        PawRAdjustmentEnd = trials_NI(j).LabelledEvents.PawRAdjustmentEnd;
        BiteBoutStart = trials_NI(j).LabelledEvents.BiteBoutStart;
        bite_timestamps = trials_NI(j).bite_timestamps;
        switch paradigm
            case 1
                temp = [MouthOpen; TongueOut; FoodinMouth; SitStart; PawLReachStart; PawRReachStart];
                if ~any((temp >= 4) & (temp <= 8)) % check if mouse drops the food between 4-8 s
                    nbite = numel(bite_timestamps(bite_timestamps > 4 & bite_timestamps <= 8));
                    if nbite == 0
                        nbite = 1;
                    end
                    nbite_NI(i) = nbite_NI(i)+nbite;
                    BiteBoutStart_NI = [BiteBoutStart_NI sum(BiteBoutStart > 4 & BiteBoutStart <= 8)];
                    Noverlap = 0;
                    for k = 1:numel(PawLAdjustmentStart)
                        if PawLAdjustmentStart(k) <= 4 && PawLAdjustmentEnd(k) > 4
                            Noverlap = Noverlap+sum(PawRAdjustmentStart > 4 & PawRAdjustmentStart < PawLAdjustmentEnd(k));
                        elseif PawLAdjustmentStart(k) > 4 && PawLAdjustmentEnd(k) <= 8
                            Noverlap = Noverlap+sum(PawRAdjustmentStart < PawLAdjustmentStart(k) & PawRAdjustmentEnd > PawLAdjustmentStart(k))+...
                                sum(PawRAdjustmentStart > PawLAdjustmentStart(k) & PawRAdjustmentStart < PawLAdjustmentEnd(k))+sum(PawRAdjustmentStart == PawLAdjustmentStart(k));
                        elseif PawLAdjustmentStart(k) <= 8 && PawLAdjustmentEnd(k) > 8
                            Noverlap = Noverlap+sum(PawRAdjustmentStart < PawLAdjustmentStart(k) & PawRAdjustmentEnd > PawLAdjustmentStart(k))+...
                                sum(PawRAdjustmentStart > PawLAdjustmentStart(k) & PawRAdjustmentStart <= 8)+sum(PawRAdjustmentStart == PawLAdjustmentStart(k));
                        end
                    end
                    PawLRAdjustment_NI = [PawLRAdjustment_NI Noverlap/nbite];
                    PawLAdjustment_NI = [PawLAdjustment_NI (sum(PawLAdjustmentStart > 4 & PawLAdjustmentStart <= 8)-Noverlap)/nbite]; % pure left
                    PawRAdjustment_NI = [PawRAdjustment_NI (sum(PawRAdjustmentStart > 4 & PawRAdjustmentStart <= 8)-Noverlap)/nbite]; % pure right
                    
                    alladj_NI(i) = alladj_NI(i)+sum(PawLAdjustmentStart > 4 & PawLAdjustmentStart <= 8)+sum(PawRAdjustmentStart > 4 & PawRAdjustmentStart <= 8);
                    biadj_NI(i) = biadj_NI(i)+Noverlap;
                end
            case 2
                if ~isempty(MouthOpen) && ~isempty(MouthRetrievalStart)
                    MouthMovement_NI = [MouthMovement_NI sum(MouthOpen < min([MouthRetrievalStart; 4]))/min([MouthRetrievalStart; 4])]; % this is event/second
                end
                if ~isempty(RetrievalStart)
                    RetrievalStartT_NI = [RetrievalStartT_NI; RetrievalStart(1)];
                end
                RetrievalTimes = numel(RetrievalStart);
                RetrievalEnd = nan(1, RetrievalTimes);
                for k = 1:RetrievalTimes
                    if RetrievalTimes == numel(SitStart)
                        if k == 1
                            MouthRetrievalStart_NI = [MouthRetrievalStart_NI sum(MouthRetrievalStart < SitEnd(k))];
                        else
                            MouthRetrievalStart_NI = [MouthRetrievalStart_NI sum(MouthRetrievalStart < SitEnd(k) & MouthRetrievalStart > SitEnd(k-1))];
                        end
                    end
                    if ~isempty(MouthRetrievalStart)
                        temp = MouthRetrievalEnd(MouthRetrievalStart == RetrievalStart(k));
                        if ~isempty(temp)
                            RetrievalEnd(k) = temp;
                        end
                        if ~isnan(RetrievalEnd(k))
                            if RetrievalTimes == numel(SitStart)
                                SitDelayT_NI = [SitDelayT_NI; SitStart(k)-RetrievalEnd(k)];
                            end
                            if RetrievalTimes == numel(PawLReachStart)
                                PawLReachDelayT_NI = [PawLReachDelayT_NI; PawLReachStart(k)-RetrievalEnd(k)];
                            end
                            if RetrievalTimes == numel(PawRReachStart)
                                PawRReachDelayT_NI = [PawRReachDelayT_NI; PawRReachStart(k)-RetrievalEnd(k)];
                            end
                        end
                    end
                end
                if ~isempty(RetrievalEnd)
                    if ~isnan(RetrievalEnd(1))
                        RetrievalEndT_NI = [RetrievalEndT_NI; RetrievalEnd(1)];
                    end
                end
                if ~isempty(SitStart)
                    SitT_NI = [SitT_NI; SitEnd-SitStart];
                    SitEndT1_NI = [SitEndT1_NI; SitEnd(1)];
                end
                if ~isempty(PawLReachStart)
                    PawLReachT_NI = [PawLReachT_NI; PawLReachEnd-PawLReachStart];
                end
                if ~isempty(PawRReachStart)
                    PawRReachT_NI = [PawRReachT_NI; PawRReachEnd-PawRReachStart];
                end
                if ~isempty(SitStart)
                    BiteDelayT_NI = [BiteDelayT_NI; bite_timestamps(1)-SitEnd(end)];
                end
                if ~isempty(SitStart)
                    PawLAdjustmentEnd(PawLAdjustmentStart < SitEnd(end) | PawLAdjustmentStart > bite_timestamps(1)) = [];
                    PawLAdjustmentStart(PawLAdjustmentStart < SitEnd(end) | PawLAdjustmentStart > bite_timestamps(1)) = [];
                    PawRAdjustmentEnd(PawRAdjustmentStart < SitEnd(end) | PawRAdjustmentStart > bite_timestamps(1)) = [];
                    PawRAdjustmentStart(PawRAdjustmentStart < SitEnd(end) | PawRAdjustmentStart > bite_timestamps(1)) = [];
                    Noverlap = 0;
                    for k = 1:numel(PawLAdjustmentStart)
                        Noverlap = Noverlap+sum(PawRAdjustmentStart < PawLAdjustmentStart(k) & PawRAdjustmentEnd > PawLAdjustmentStart(k))+...
                                sum(PawRAdjustmentStart > PawLAdjustmentStart(k) & PawRAdjustmentStart < PawLAdjustmentEnd(k))+sum(PawRAdjustmentStart == PawLAdjustmentStart(k));
                    end
                    PawLRAdjustment_NI = [PawLRAdjustment_NI Noverlap];
                    PawLAdjustment_NI = [PawLAdjustment_NI (numel(PawLAdjustmentStart)-Noverlap)];
                    PawRAdjustment_NI = [PawRAdjustment_NI (numel(PawRAdjustmentStart)-Noverlap)];
                end
        end
    end
    APawLAdjustment_NI(i) = mean(PawLAdjustment_NI);
    APawRAdjustment_NI(i) = mean(PawRAdjustment_NI);
    APawLRAdjustment_NI(i) = mean(PawLRAdjustment_NI);
    
    AllPawLAdjustment_NI = [AllPawLAdjustment_NI PawLAdjustment_NI];
    AllPawRAdjustment_NI = [AllPawRAdjustment_NI PawRAdjustment_NI];
    AllPawLRAdjustment_NI = [AllPawLRAdjustment_NI PawLRAdjustment_NI];
    switch paradigm
        case 1
            ABiteBoutStart_NI(i) = mean(BiteBoutStart_NI);
        case 2
            AMouthMovement_NI(i) = mean(MouthMovement_NI);
            ARetrievalStartT_NI(i) = mean(RetrievalStartT_NI);
            ARetrievalEndT_NI(i) = mean(RetrievalEndT_NI);
            AMouthRetrievalStart_NI(i) = mean(MouthRetrievalStart_NI);
            ASitDelayT_NI(i) = mean(SitDelayT_NI);
            APawLReachDelayT_NI(i) = mean(PawLReachDelayT_NI);
            APawRReachDelayT_NI(i) = mean(PawRReachDelayT_NI);
            ASitT_NI(i) = mean(SitT_NI);
            APawLReachT_NI(i) = mean(PawLReachT_NI);
            APawRReachT_NI(i) = mean(PawRReachT_NI);
            ABiteDelayT_NI(i) = mean(BiteDelayT_NI);
            ASitEndT1_NI(i) = mean(SitEndT1_NI);
    end
    
    for j = 1:numel(trials_I)
        MouthOpen = trials_I(j).LabelledEvents.MouthOpen;
        TongueOut = trials_I(j).LabelledEvents.TongueOut;
        FoodinMouth = trials_I(j).LabelledEvents.FoodinMouth;
        SitStart = trials_I(j).LabelledEvents.SitStart;
        SitEnd = trials_I(j).LabelledEvents.SitEnd;
        PawLReachStart = trials_I(j).LabelledEvents.PawLReachStart;
        PawLReachEnd = trials_I(j).LabelledEvents.PawLReachEnd;
        PawRReachStart = trials_I(j).LabelledEvents.PawRReachStart;
        PawRReachEnd = trials_I(j).LabelledEvents.PawRReachEnd;
        MouthRetrievalStart = trials_I(j).LabelledEvents.MouthRetrievalStart;
        MouthRetrievalEnd = trials_I(j).LabelledEvents.MouthRetrievalEnd;
        RetrievalStart = trials_I(j).LabelledEvents.RetrievalStart;
        PawLAdjustmentStart = trials_I(j).LabelledEvents.PawLAdjustmentStart;
        PawLAdjustmentEnd = trials_I(j).LabelledEvents.PawLAdjustmentEnd;
        PawRAdjustmentStart = trials_I(j).LabelledEvents.PawRAdjustmentStart;
        PawRAdjustmentEnd = trials_I(j).LabelledEvents.PawRAdjustmentEnd;
        BiteBoutStart = trials_I(j).LabelledEvents.BiteBoutStart;
        bite_timestamps = trials_I(j).bite_timestamps;
        
        switch paradigm
            case 1
                temp = [MouthOpen; TongueOut; FoodinMouth; SitStart; PawLReachStart; PawRReachStart];
                if ~any((temp >= 4) & (temp <= 8)) % check if mouse drops the food between 4-8 s
                    nbite = numel(bite_timestamps(bite_timestamps > 4 & bite_timestamps <= 8));
                    if nbite == 0
                        nbite = 1;
                    end
                    nbite_I(i) = nbite_I(i)+nbite;
                    BiteBoutStart_I = [BiteBoutStart_I sum(BiteBoutStart > 4 & BiteBoutStart <= 8)];
                    Noverlap = 0;
                    for k = 1:numel(PawLAdjustmentStart)
                        if PawLAdjustmentStart(k) <= 4 && PawLAdjustmentEnd(k) > 4
                            Noverlap = Noverlap+sum(PawRAdjustmentStart > 4 & PawRAdjustmentStart < PawLAdjustmentEnd(k));
                        elseif PawLAdjustmentStart(k) > 4 && PawLAdjustmentEnd(k) <= 8
                            Noverlap = Noverlap+sum(PawRAdjustmentStart < PawLAdjustmentStart(k) & PawRAdjustmentEnd > PawLAdjustmentStart(k))+...
                                sum(PawRAdjustmentStart > PawLAdjustmentStart(k) & PawRAdjustmentStart < PawLAdjustmentEnd(k))+sum(PawRAdjustmentStart == PawLAdjustmentStart(k));
                        elseif PawLAdjustmentStart(k) > 4 && PawLAdjustmentStart(k) <= 8 && PawLAdjustmentEnd(k) > 8
                            Noverlap = Noverlap+sum(PawRAdjustmentStart < PawLAdjustmentStart(k) & PawRAdjustmentEnd > PawLAdjustmentStart(k))+...
                                sum(PawRAdjustmentStart > PawLAdjustmentStart(k) & PawRAdjustmentStart <= 8)+sum(PawRAdjustmentStart == PawLAdjustmentStart(k));
                        end
                    end
                    PawLRAdjustment_I = [PawLRAdjustment_I Noverlap/nbite];
                    PawLAdjustment_I = [PawLAdjustment_I (sum(PawLAdjustmentStart > 4 & PawLAdjustmentStart <= 8)-Noverlap)/nbite];
                    PawRAdjustment_I = [PawRAdjustment_I (sum(PawRAdjustmentStart > 4 & PawRAdjustmentStart <= 8)-Noverlap)/nbite];
                    
                    alladj_I(i) = alladj_I(i)+sum(PawLAdjustmentStart > 4 & PawLAdjustmentStart <= 8)+sum(PawRAdjustmentStart > 4 & PawRAdjustmentStart <= 8);
                    biadj_I(i) = biadj_I(i)+Noverlap;
                end
            case 2
                if ~isempty(MouthOpen) && ~isempty(MouthRetrievalStart)
                    MouthMovement_I = [MouthMovement_I sum(MouthOpen < min([MouthRetrievalStart; 4]))/min([MouthRetrievalStart; 4])]; % this is event/second
                end
                if ~isempty(RetrievalStart)
                    RetrievalStartT_I = [RetrievalStartT_I; RetrievalStart(1)];
                end
                RetrievalTimes = numel(RetrievalStart);
                RetrievalEnd = nan(1, RetrievalTimes);
                for k = 1:RetrievalTimes
                    if RetrievalTimes == numel(SitStart)
                        if k == 1
                            temp = MouthRetrievalStart(MouthRetrievalStart < SitEnd(k));
                        else
                            temp = MouthRetrievalStart(MouthRetrievalStart < SitEnd(k) & MouthRetrievalStart > SitEnd(k-1));
                        end
                        if any(temp < 4)
                            MouthRetrievalStart_I = [MouthRetrievalStart_I numel(temp)];
                        end
                    end
                    if ~isempty(MouthRetrievalStart)
                        temp = MouthRetrievalEnd(MouthRetrievalStart == RetrievalStart(k));
                        if ~isempty(temp)
                            RetrievalEnd(k) = temp;
                        end
                        if ~isnan(RetrievalEnd(k))
                            if RetrievalEnd(k) < 4
                                if RetrievalTimes == numel(SitStart)
                                    SitDelayT_I = [SitDelayT_I; SitStart(k)-RetrievalEnd(k)];
                                end
                                if RetrievalTimes == numel(PawLReachStart)
                                    PawLReachDelayT_I = [PawLReachDelayT_I; PawLReachStart(k)-RetrievalEnd(k)];
                                end
                                if RetrievalTimes == numel(PawRReachStart)
                                    PawRReachDelayT_I = [PawRReachDelayT_I; PawRReachStart(k)-RetrievalEnd(k)];
                                end
                            end
                        end
                    end
                    if ~isempty(RetrievalEnd)
                        if ~isnan(RetrievalEnd(1))
                            RetrievalEndT_I = [RetrievalEndT_I; RetrievalEnd(1)];
                        end
                    end
                    if ~isempty(SitStart)
                        try
                            if SitStart(k) < 4
                                SitT_I = [SitT_I; SitEnd(k)-SitStart(k)];
                            end
                        end
                        SitEndT1_I = [SitEndT1_I; SitEnd(1)];
                    end
                    if ~isempty(PawLReachStart)
                        try
                            if PawLReachStart(k) < 4
                                PawLReachT_I = [PawLReachT_I; PawLReachEnd(k)-PawLReachStart(k)];
                            end
                        end
                    end
                    if ~isempty(PawRReachStart)
                        try
                            if PawRReachStart(k) < 4
                                PawRReachT_I = [PawRReachT_I; PawRReachEnd(k)-PawRReachStart(k)];
                            end
                        end
                    end
                end
                if ~isempty(SitStart)
                    if SitEnd(end) < 4
                        BiteDelayT_I = [BiteDelayT_I; bite_timestamps(1)-SitEnd(end)];
                    end
                end
                if ~isempty(SitStart)
                    if SitEnd(end) < 4
                        PawLAdjustmentEnd(PawLAdjustmentStart < SitEnd(end) | PawLAdjustmentStart > bite_timestamps(1)) = [];
                        PawLAdjustmentStart(PawLAdjustmentStart < SitEnd(end) | PawLAdjustmentStart > bite_timestamps(1)) = [];
                        PawRAdjustmentEnd(PawRAdjustmentStart < SitEnd(end) | PawRAdjustmentStart > bite_timestamps(1)) = [];
                        PawRAdjustmentStart(PawRAdjustmentStart < SitEnd(end) | PawRAdjustmentStart > bite_timestamps(1)) = [];
                        Noverlap = 0;
                        for k = 1:numel(PawLAdjustmentStart)
                            Noverlap = Noverlap+sum(PawRAdjustmentStart < PawLAdjustmentStart(k) & PawRAdjustmentEnd > PawLAdjustmentStart(k))+...
                                sum(PawRAdjustmentStart > PawLAdjustmentStart(k) & PawRAdjustmentStart < PawLAdjustmentEnd(k))+sum(PawRAdjustmentStart == PawLAdjustmentStart(k));
                        end
                        PawLRAdjustment_I = [PawLRAdjustment_I Noverlap];
                        PawLAdjustment_I = [PawLAdjustment_I (numel(PawLAdjustmentStart)-Noverlap)];
                        PawRAdjustment_I = [PawRAdjustment_I (numel(PawRAdjustmentStart)-Noverlap)];
                    end
                end
        end
    end
    APawLAdjustment_I(i) = mean(PawLAdjustment_I);
    APawRAdjustment_I(i) = mean(PawRAdjustment_I);
    APawLRAdjustment_I(i) = mean(PawLRAdjustment_I);
    
    AllPawLAdjustment_I = [AllPawLAdjustment_I PawLAdjustment_I];
    AllPawRAdjustment_I = [AllPawRAdjustment_I PawRAdjustment_I];
    AllPawLRAdjustment_I = [AllPawLRAdjustment_I PawLRAdjustment_I];
    switch paradigm
        case 1
            ABiteBoutStart_I(i) = mean(BiteBoutStart_I);
        case 2
            AMouthMovement_I(i) = mean(MouthMovement_I);
            ARetrievalStartT_I(i) = mean(RetrievalStartT_I);
            ARetrievalEndT_I(i) = mean(RetrievalEndT_I);
            AMouthRetrievalStart_I(i) = mean(MouthRetrievalStart_I);
            ASitDelayT_I(i) = mean(SitDelayT_I);
            APawLReachDelayT_I(i) = mean(PawLReachDelayT_I);
            APawRReachDelayT_I(i) = mean(PawRReachDelayT_I);
            ASitT_I(i) = mean(SitT_I);
            APawLReachT_I(i) = mean(PawLReachT_I);
            APawRReachT_I(i) = mean(PawRReachT_I);
            ABiteDelayT_I(i) = mean(BiteDelayT_I);
            ASitEndT1_I(i) = mean(SitEndT1_I);
    end
end
switch paradigm
    case 1
        ABiteBoutStart_I = ABiteBoutStart_I/4;
        ABiteBoutStart_NI = ABiteBoutStart_NI/4;
        
        figure;
        paired_plot(APawLAdjustment_NI+APawRAdjustment_NI+APawLRAdjustment_NI, APawLAdjustment_I+APawRAdjustment_I+APawLRAdjustment_I, 'Adjustments/bite', filename);
        figure;
        paired_plot(APawLRAdjustment_NI, APawLRAdjustment_I, 'Bimanual Adjustments/bite', filename);
        figure;
        paired_plot(APawLRAdjustment_NI./(APawLAdjustment_NI+APawRAdjustment_NI+APawLRAdjustment_NI),...
            APawLRAdjustment_I./(APawLAdjustment_I+APawRAdjustment_I+APawLRAdjustment_I), 'Bimanual/All', filename);
        figure;
        paired_plot(ABiteBoutStart_NI, ABiteBoutStart_I, 'Bouts/s', filename);
        
%         figure;
%         paired_plot((alladj_NI-biadj_NI)./nbite_NI, (alladj_I-biadj_I)./nbite_I, 'adjustment/bite', filename);
%         figure;
%         paired_plot(biadj_NI./nbite_NI, biadj_I./nbite_I, 'bimanual adjustment/bite', filename);
%         figure;
%         paired_plot(biadj_NI./(alladj_NI-biadj_NI), biadj_I./(alladj_I-biadj_I), 'bimanual/all', filename);
        
        compare2distributions(AllPawLAdjustment_NI+AllPawRAdjustment_NI+AllPawLRAdjustment_NI, AllPawLAdjustment_I+AllPawRAdjustment_I+AllPawLRAdjustment_I, 'Adjustments/bite', '');
        compare2distributions(AllPawLRAdjustment_NI, AllPawLRAdjustment_I, 'Bimanual Adjustments/bite', '');
        
        violin_boxplot(AllPawLAdjustment_NI+AllPawRAdjustment_NI+AllPawLRAdjustment_NI, AllPawLAdjustment_I+AllPawRAdjustment_I+AllPawLRAdjustment_I, 'Adjustments/bite');
        violin_boxplot(AllPawLRAdjustment_NI, AllPawLRAdjustment_I, 'Bimanual Adjustments/bite');
    case 2
        figure;
        paired_plot(AMouthMovement_NI, AMouthMovement_I, 'Jaw Movement/s', filename);
        figure;
        paired_plot(ARetrievalStartT_NI, ARetrievalStartT_I, 'Retrieval Success Time (s)', filename); % retrieval with both mouth and hand
        figure;
        paired_plot(ARetrievalEndT_NI, ARetrievalEndT_I, 'Retrieval End Time (s)', filename); % only consider retrieval with mouth
        figure;
        paired_plot(AMouthRetrievalStart_NI, AMouthRetrievalStart_I, '# of Jaw Retrieval', filename);
        figure;
        paired_plot(ASitDelayT_NI*10^3, ASitDelayT_I*10^3, 'Delay to Sit (ms)', filename, 'dot');
        figure;
        paired_plot(APawLReachDelayT_NI*10^3, APawLReachDelayT_I*10^3, 'Delay to PawLReach (ms)', filename, 'dot');
        figure;
        paired_plot(APawRReachDelayT_NI*10^3, APawRReachDelayT_I*10^3, 'Delay to PawRReach (ms)', filename, 'dot');
        figure;
        paired_plot(ASitT_NI*10^3, ASitT_I*10^3, 'Sit (ms)', filename);
        figure;
        paired_plot(APawLReachT_NI*10^3, APawLReachT_I*10^3, 'PawLReach (ms)', filename);
        figure;
        paired_plot(APawRReachT_NI*10^3, APawRReachT_I*10^3, 'PawRReach (ms)', filename);
        figure;
        paired_plot(ABiteDelayT_NI, ABiteDelayT_I, 'Delay to 1st Bite (s)', filename);
        figure;
        paired_plot(APawLAdjustment_NI+APawRAdjustment_NI+APawLRAdjustment_NI, APawLAdjustment_I+APawRAdjustment_I+APawLRAdjustment_I, '# of Adjustments', filename);
        figure;
        paired_plot(APawLRAdjustment_NI, APawLRAdjustment_I, '# of Bimanual Adjustments', filename);
        figure;
        paired_plot(ASitEndT1_NI, ASitEndT1_I, 'Time of 1st Sit End (s)', filename);
end

%% for hand adjustment duration and delay of first hand adjustment to hand withdraw
[filename, pathname] = uigetfile('*.mat', 'Pick all of the individual data set', 'MultiSelect', 'on');
if isequal(filename, 0)
    return;
end
clc;
N = numel(filename);
PawAdj_I = cell(1, N);
PawAdj_NI = cell(1, N);
PawAdj2withdraw_I = cell(1, N);
PawAdj2withdraw_NI = cell(1, N);

for i = 1:N
    PawAdj_duration = [];
    PawAdj2withdraw = [];
    
    load([pathname filename{i}]);
    for j = 1:numel(trials_NI)
        MouthOpen = trials_NI(j).LabelledEvents.MouthOpen;
        TongueOut = trials_NI(j).LabelledEvents.TongueOut;
        FoodinMouth = trials_NI(j).LabelledEvents.FoodinMouth;
        SitStart = trials_NI(j).LabelledEvents.SitStart;
        SitEnd = trials_NI(j).LabelledEvents.SitEnd;
        PawLReachStart = trials_NI(j).LabelledEvents.PawLReachStart;
        PawLReachEnd = trials_NI(j).LabelledEvents.PawLReachEnd;
        PawRReachStart = trials_NI(j).LabelledEvents.PawRReachStart;
        PawRReachEnd = trials_NI(j).LabelledEvents.PawRReachEnd;
        MouthRetrievalStart = trials_NI(j).LabelledEvents.MouthRetrievalStart;
        MouthRetrievalEnd = trials_NI(j).LabelledEvents.MouthRetrievalEnd;
        RetrievalStart = trials_NI(j).LabelledEvents.RetrievalStart;
        
        PawLAdjustmentStart = trials_NI(j).LabelledEvents.PawLAdjustmentStart;
        PawLAdjustmentEnd = trials_NI(j).LabelledEvents.PawLAdjustmentEnd;
        PawRAdjustmentStart = trials_NI(j).LabelledEvents.PawRAdjustmentStart;
        PawRAdjustmentEnd = trials_NI(j).LabelledEvents.PawRAdjustmentEnd;
        BiteBoutStart = trials_NI(j).LabelledEvents.BiteBoutStart;
        BiteBoutStart = sort(BiteBoutStart);
        
        temp = [MouthOpen; TongueOut; FoodinMouth; SitStart; PawLReachStart; PawRReachStart];
        if ~any((temp >= 4) & (temp <= 8)) % check if mouse drops the food between 4-8 s
            PawLAdjustmentStart_valid = PawLAdjustmentStart(PawLAdjustmentStart >= 4 & PawLAdjustmentStart < 8);
            PawLAdjustmentEnd_valid = PawLAdjustmentEnd(PawLAdjustmentStart >= 4 & PawLAdjustmentStart < 8);
            PawRAdjustmentStart_valid = PawRAdjustmentStart(PawRAdjustmentStart >= 4 & PawRAdjustmentStart < 8);
            PawRAdjustmentEnd_valid = PawRAdjustmentEnd(PawRAdjustmentStart >= 4 & PawRAdjustmentStart < 8);
            PawAdj_duration = [PawAdj_duration; PawLAdjustmentEnd_valid-PawLAdjustmentStart_valid; PawRAdjustmentEnd_valid-PawRAdjustmentStart_valid];
            
            BiteBoutStart_valid = BiteBoutStart(BiteBoutStart >= 4 & BiteBoutStart < 8);
            BiteBoutStart = BiteBoutStart(BiteBoutStart >= 8);
            BiteBoutStart_valid = [BiteBoutStart_valid; BiteBoutStart(1)];
            if numel(BiteBoutStart_valid) > 1
                for k = 1:numel(BiteBoutStart_valid)-1
                    PawAdj2withdraw = [PawAdj2withdraw; min([PawLAdjustmentStart(PawLAdjustmentStart > BiteBoutStart_valid(k) & PawLAdjustmentStart < BiteBoutStart_valid(k+1));...
                        PawRAdjustmentStart(PawRAdjustmentStart > BiteBoutStart_valid(k) & PawRAdjustmentStart < BiteBoutStart_valid(k+1))]-BiteBoutStart_valid(k))];
                end
            end
        end
    end
    PawAdj_NI{i} = PawAdj_duration;
    PawAdj2withdraw_NI{i} = PawAdj2withdraw;
    
    PawAdj_duration = [];
    PawAdj2withdraw = [];
    for j = 1:numel(trials_I)
        MouthOpen = trials_I(j).LabelledEvents.MouthOpen;
        TongueOut = trials_I(j).LabelledEvents.TongueOut;
        FoodinMouth = trials_I(j).LabelledEvents.FoodinMouth;
        SitStart = trials_I(j).LabelledEvents.SitStart;
        SitEnd = trials_I(j).LabelledEvents.SitEnd;
        PawLReachStart = trials_I(j).LabelledEvents.PawLReachStart;
        PawLReachEnd = trials_I(j).LabelledEvents.PawLReachEnd;
        PawRReachStart = trials_I(j).LabelledEvents.PawRReachStart;
        PawRReachEnd = trials_I(j).LabelledEvents.PawRReachEnd;
        MouthRetrievalStart = trials_I(j).LabelledEvents.MouthRetrievalStart;
        MouthRetrievalEnd = trials_I(j).LabelledEvents.MouthRetrievalEnd;
        RetrievalStart = trials_I(j).LabelledEvents.RetrievalStart;
        
        PawLAdjustmentStart = trials_I(j).LabelledEvents.PawLAdjustmentStart;
        PawLAdjustmentEnd = trials_I(j).LabelledEvents.PawLAdjustmentEnd;
        PawRAdjustmentStart = trials_I(j).LabelledEvents.PawRAdjustmentStart;
        PawRAdjustmentEnd = trials_I(j).LabelledEvents.PawRAdjustmentEnd;
        BiteBoutStart = trials_I(j).LabelledEvents.BiteBoutStart;
        
        temp = [MouthOpen; TongueOut; FoodinMouth; SitStart; PawLReachStart; PawRReachStart];
        if ~any((temp >= 4) & (temp <= 8)) % check if mouse drops the food between 4-8 s
            PawLAdjustmentStart_valid = PawLAdjustmentStart(PawLAdjustmentStart >= 4 & PawLAdjustmentStart < 8);
            PawLAdjustmentEnd_valid = PawLAdjustmentEnd(PawLAdjustmentStart >= 4 & PawLAdjustmentStart < 8);
            PawRAdjustmentStart_valid = PawRAdjustmentStart(PawRAdjustmentStart >= 4 & PawRAdjustmentStart < 8);
            PawRAdjustmentEnd_valid = PawRAdjustmentEnd(PawRAdjustmentStart >= 4 & PawRAdjustmentStart < 8);
            PawAdj_duration = [PawAdj_duration; PawLAdjustmentEnd_valid-PawLAdjustmentStart_valid; PawRAdjustmentEnd_valid-PawRAdjustmentStart_valid];
            
            BiteBoutStart_valid = BiteBoutStart(BiteBoutStart >= 4 & BiteBoutStart < 8);
            BiteBoutStart = BiteBoutStart(BiteBoutStart >= 8);
            BiteBoutStart_valid = [BiteBoutStart_valid; BiteBoutStart(1)];
            if numel(BiteBoutStart_valid) > 1
                for k = 1:numel(BiteBoutStart_valid)-1
                    PawAdj2withdraw = [PawAdj2withdraw; min([PawLAdjustmentStart(PawLAdjustmentStart > BiteBoutStart_valid(k) & PawLAdjustmentStart < BiteBoutStart_valid(k+1));...
                        PawRAdjustmentStart(PawRAdjustmentStart > BiteBoutStart_valid(k) & PawRAdjustmentStart < BiteBoutStart_valid(k+1))]-BiteBoutStart_valid(k))];
                end
            end
        end
    end
    PawAdj_I{i} = PawAdj_duration;
    PawAdj2withdraw_I{i} = PawAdj2withdraw;
end

PawAdj_NI_avg = nan(1, N);
PawAdj_I_avg = nan(1, N);
PawAdj_NI_all = [];
PawAdj_I_all = [];
PawAdj2withdraw_NI_avg = nan(1, N);
PawAdj2withdraw_I_avg = nan(1, N);
PawAdj2withdraw_NI_all = [];
PawAdj2withdraw_I_all = [];
for i = 1:N
    PawAdj_NI_avg(i) = mean(PawAdj_NI{i});
    PawAdj_NI_all = [PawAdj_NI_all; PawAdj_NI{i}];
    PawAdj_I_avg(i) = mean(PawAdj_I{i});
    PawAdj_I_all = [PawAdj_I_all; PawAdj_I{i}];
    PawAdj2withdraw_NI_avg(i) = mean(PawAdj2withdraw_NI{i});
    PawAdj2withdraw_NI_all = [PawAdj2withdraw_NI_all; PawAdj2withdraw_NI{i}];
    PawAdj2withdraw_I_avg(i) = mean(PawAdj2withdraw_I{i});
    PawAdj2withdraw_I_all = [PawAdj2withdraw_I_all; PawAdj2withdraw_I{i}];
end

figure;
paired_plot(PawAdj_NI_avg, PawAdj_I_avg, 'Adjustment duration (s)', filename);
figure;
paired_plot(PawAdj2withdraw_NI_avg, PawAdj2withdraw_I_avg, '1st Adjustment to withdraw (s)', filename);

compare2distributions(PawAdj_NI_all, PawAdj_I_all, 'Adjustment duration (s)', '');
compare2distributions(PawAdj2withdraw_NI_all, PawAdj2withdraw_I_all, '1st Adjustment to withdraw (s)', '');

violin_boxplot(PawAdj_NI_all, PawAdj_I_all, 'Adjustment duration (s)');
violin_boxplot(PawAdj2withdraw_NI_all, PawAdj2withdraw_I_all, '1st Adjustment to withdraw (s)');

resampling = 10000;
figure;
permutation_test(PawAdj_NI, PawAdj_I, 'both', resampling, 'Adjustment duration (s)');
figure;
permutation_test(PawAdj2withdraw_NI, PawAdj2withdraw_I, 'both', resampling, '1st Adjustment to withdraw (s)');
figure;
bootstrapping_test(PawAdj_NI, PawAdj_I, 'both', resampling, 'Adjustment duration (s)');
figure;
bootstrapping_test(PawAdj2withdraw_NI, PawAdj2withdraw_I, 'both', resampling, '1st Adjustment to withdraw (s)');