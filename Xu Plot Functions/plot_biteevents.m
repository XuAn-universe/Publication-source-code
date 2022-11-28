function ybite = plot_biteevents(timestamps, amplitudes)
legend(gca, 'off');
hl = findobj(gca, 'Type', 'line');
yl = get(hl, 'YData');
tl = get(hl, 'XData');
yrange = range(yl);
ymax = max(yl);
ymin = min(yl);
scalefactor = 0.1;
ybite = [];
for j = 1:numel(timestamps)
    timestamp = timestamps(j);
    if timestamp >= tl(1) && timestamp <= tl(end)
        [~, id] = min(abs(tl-timestamp));
        line([timestamp timestamp], [yl(id)+yrange*scalefactor*amplitudes(j) yl(id)-yrange*scalefactor*amplitudes(j)],...
            'Color', [1 0 0], 'LineStyle', '-', 'LineWidth', 1);
        ymax = max([ymax yl(id)+yrange*scalefactor*amplitudes(j)]);
        ymin = min([ymin yl(id)-yrange*scalefactor*amplitudes(j)]);
        ybite = [ybite yl(id)];
    end
end
% ylim([ymin ymax]);
end