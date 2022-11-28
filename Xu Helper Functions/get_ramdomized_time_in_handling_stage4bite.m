function bite_timestampsRnd = get_ramdomized_time_in_handling_stage4bite(events, bite_timestamps)
RetrievalStart = events.RetrievalStart;
BiteBoutStart = events.BiteBoutStart;
BiteBoutStart = [RetrievalStart; BiteBoutStart; inf];
BiteBoutStart = sort(BiteBoutStart);
bite_timestampsRnd = nan(size(bite_timestamps));
for i = 1:numel(bite_timestamps)
    BoutStart = BiteBoutStart(find(BiteBoutStart < bite_timestamps(i), 1, 'last'));
    if ~(isempty(events.FeedingEnd) && BoutStart == BiteBoutStart(end-1))
        BoutEnd = BiteBoutStart(find(BiteBoutStart > bite_timestamps(i), 1, 'first'));
        bitetime = bite_timestamps(bite_timestamps > BoutStart & bite_timestamps < BoutEnd);
        chewstart = bitetime(end);
        bite_timestampsRnd(i) = BoutStart+(chewstart-BoutStart)*rand(1, 1);
    end
end
bite_timestampsRnd(isnan(bite_timestampsRnd)) = [];