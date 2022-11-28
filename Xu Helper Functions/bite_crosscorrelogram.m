function bite_crosscorrelogram(bite_events, bite_events_I)
% By Xu An Nov. 2019 @Cold Spring Harbor Laboratory
% xan@cshl.edu
tmax = 10;
resolution = 0.1;
trials = numel(bite_events);
interval_all = cell(1, trials);
intervals = [];
for i = 1:trials
    t = bite_events(i).bite_timestamps;
    if ~isempty(t)
        interval = t'*ones(1, numel(t))-ones(numel(t), 1)*t;
        interval_all{i} = interval(:);
        intervals = [intervals; interval(:)];
    end
end
figure;
histogram(intervals(intervals >= -tmax & intervals <= tmax), linspace(-tmax, tmax, 2*tmax/resolution), 'FaceColor', [0 0 0], 'EdgeColor', [1 1 1], 'FaceAlpha', 1, 'Normalization', 'probability');
hold on;
if nargin == 2
    trials = numel(bite_events_I);
    interval_all = cell(1, trials);
    intervals = [];
    for i = 1:trials
        t = bite_events_I(i).bite_timestamps;
        if ~isempty(t)
            interval = t'*ones(1, numel(t))-ones(numel(t), 1)*t;
            interval_all{i} = interval(:);
            intervals = [intervals; interval(:)];
        end
    end
    histogram(intervals(intervals >= -tmax & intervals <= tmax), linspace(-tmax, tmax, 2*tmax/resolution), 'FaceColor', 'none', 'EdgeColor', [0 1 0], 'FaceAlpha', 1, 'Normalization', 'probability');
end
set(gca, 'TickLength', [0 0], 'FontSize', 12);
box off;
xlabel('Bite Interval (s)');
ylabel('Probability');
title('Bite Crosscorrelogram')