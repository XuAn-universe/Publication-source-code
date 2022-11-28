function result = get_ramdomized_time_in_handling_stage(events, bite_timestamps)
PawLAdjustmentStart = events.PawLAdjustmentStart;
PawRAdjustmentStart = events.PawRAdjustmentStart;
RetrievalStart = events.RetrievalStart;
BiteBoutStart = events.BiteBoutStart;
BiteBoutStart = [RetrievalStart; BiteBoutStart; inf];
BiteBoutStart = sort(BiteBoutStart);
nobiteID = [];
PawLAdjustmentRnd = [];
for i = 1:numel(PawLAdjustmentStart)
    BoutStart = BiteBoutStart(find(BiteBoutStart < PawLAdjustmentStart(i), 1, 'last'));
    BoutEnd = BiteBoutStart(find(BiteBoutStart > PawLAdjustmentStart(i), 1, 'first'));
    bitetime = bite_timestamps(bite_timestamps > BoutStart & bite_timestamps < BoutEnd);
    if isempty(bitetime)
        nobiteID = [nobiteID i];
    else
        chewstart = bitetime(end);
        PawLAdjustmentRnd = [PawLAdjustmentRnd BoutStart+(chewstart-BoutStart)*rand(1, 1)];
    end
end
PawLAdjustmentStart(nobiteID) = [];
result.PawLAdjustmentStart = PawLAdjustmentStart;
result.PawLAdjustmentRnd = PawLAdjustmentRnd;

nobiteID = [];
PawRAdjustmentRnd = [];
for i = 1:numel(PawRAdjustmentStart)
    BoutStart = BiteBoutStart(find(BiteBoutStart < PawRAdjustmentStart(i), 1, 'last'));
    BoutEnd = BiteBoutStart(find(BiteBoutStart > PawRAdjustmentStart(i), 1, 'first'));
    bitetime = bite_timestamps(bite_timestamps > BoutStart & bite_timestamps < BoutEnd);
    if isempty(bitetime)
        nobiteID = [nobiteID i];
    else
        chewstart = bitetime(end);
        PawRAdjustmentRnd = [PawRAdjustmentRnd BoutStart+(chewstart-BoutStart)*rand(1, 1)];
    end
end
PawRAdjustmentStart(nobiteID) = [];
result.PawRAdjustmentStart = PawRAdjustmentStart;
result.PawRAdjustmentRnd = PawRAdjustmentRnd;