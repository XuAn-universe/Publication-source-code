function plot_events(hsp, i, trials, LabelledEvents)
for j = 1:numel(LabelledEvents.PawLReachStart)
    patch(hsp, [LabelledEvents.PawLReachStart(j) LabelledEvents.PawLReachStart(j) LabelledEvents.PawLReachEnd(j), LabelledEvents.PawLReachEnd(j)],...
        [(trials-i)*6+1 (trials-i)*6+2 (trials-i)*6+2 (trials-i)*6+1], [0 0 0.5], 'EdgeColor', 'none');
end
for j = 1:numel(LabelledEvents.PawLAdjustmentStart)
    patch(hsp, [LabelledEvents.PawLAdjustmentStart(j) LabelledEvents.PawLAdjustmentStart(j) LabelledEvents.PawLAdjustmentEnd(j), LabelledEvents.PawLAdjustmentEnd(j)],...
        [(trials-i)*6+1 (trials-i)*6+2 (trials-i)*6+2 (trials-i)*6+1], [0 0 1], 'EdgeColor', 'none');
end
for j = 1:numel(LabelledEvents.PawRReachStart)
    patch(hsp, [LabelledEvents.PawRReachStart(j) LabelledEvents.PawRReachStart(j) LabelledEvents.PawRReachEnd(j), LabelledEvents.PawRReachEnd(j)],...
        [(trials-i)*6+2 (trials-i)*6+3 (trials-i)*6+3 (trials-i)*6+2], [0 0.5 0], 'EdgeColor', 'none');
end
for j = 1:numel(LabelledEvents.PawRAdjustmentStart)
    patch(hsp, [LabelledEvents.PawRAdjustmentStart(j) LabelledEvents.PawRAdjustmentStart(j) LabelledEvents.PawRAdjustmentEnd(j), LabelledEvents.PawRAdjustmentEnd(j)],...
        [(trials-i)*6+2 (trials-i)*6+3 (trials-i)*6+3 (trials-i)*6+2], [0 1 0], 'EdgeColor', 'none');
end
for j = 1:numel(LabelledEvents.SitStart)
    patch(hsp, [LabelledEvents.SitStart(j) LabelledEvents.SitStart(j) LabelledEvents.SitEnd(j), LabelledEvents.SitEnd(j)],...
        [(trials-i)*6+3 (trials-i)*6+4 (trials-i)*6+4 (trials-i)*6+3], [0.4940 0.1840 0.5560], 'EdgeColor', 'none');
end
for j = 1:numel(LabelledEvents.MouthOpen)
    if any(LabelledEvents.MouthRetrievalStart == LabelledEvents.MouthOpen(j)) %%% uncomment to show all jaw movements
        patch(hsp, [LabelledEvents.MouthOpen(j) LabelledEvents.MouthOpen(j) LabelledEvents.MouthClosed(j), LabelledEvents.MouthClosed(j)],...
            [(trials-i)*6+5 (trials-i)*6+6 (trials-i)*6+6 (trials-i)*6+5], [0.6350 0.0780 0.1840], 'EdgeColor', 'none');
    end %%%
end
% uncomment to show all jaw movements
% for j = 1:numel(LabelledEvents.MouthRetrievalStart)
%     plot(hsp, (LabelledEvents.MouthRetrievalStart(j)+LabelledEvents.MouthRetrievalEnd(j))/2, (trials-i)*6+5.5, '*k', 'MarkerSize', 12);
% end

for j = 1:numel(LabelledEvents.FoodinMouth)
    line(hsp, [LabelledEvents.FoodinMouth(j) LabelledEvents.FoodinMouth(j)], [(trials-i)*6+4 (trials-i)*6+6],...
        'Color', [0.4660 0.6740 0.1880], 'LineStyle', '-', 'LineWidth', 1);
end
for j = 1:numel(LabelledEvents.TongueOut)
    patch(hsp, [LabelledEvents.TongueOut(j) LabelledEvents.TongueOut(j) LabelledEvents.TongueIn(j), LabelledEvents.TongueIn(j)],...
        [(trials-i)*6+4 (trials-i)*6+5 (trials-i)*6+5 (trials-i)*6+4], [0.9290 0.6940 0.1250], 'EdgeColor', 'none');
end
for j = 1:numel(LabelledEvents.BiteBoutStart)
    line(hsp, [LabelledEvents.BiteBoutStart(j) LabelledEvents.BiteBoutStart(j)], [(trials-i)*6 (trials-i)*6+3],...
        'Color', [0 1 1], 'LineStyle', '-', 'LineWidth', 1);
end
for j = 1:numel(LabelledEvents.FeedingEnd)
    line(hsp, [LabelledEvents.FeedingEnd(j) LabelledEvents.FeedingEnd(j)], [(trials-i)*6 (trials-i)*6+6],...
        'Color', [1 0 0], 'LineStyle', '-', 'LineWidth', 2);
end