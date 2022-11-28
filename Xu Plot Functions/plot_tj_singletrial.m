function hp = plot_tj_singletrial(t, y, time_range, laserstart, laserstop, pcolor, lcolor, ylabel_text, llabel)
scalefactor = 0.1;
y = y(t >= time_range(1) & t <= time_range(2), :);
t = t(t >= time_range(1) & t <= time_range(2));
yrange = range(y(:));
hold on;
if ~isempty(laserstart)
    for j = 1:numel(laserstart)
%         patch([laserstart(j) laserstart(j) laserstop(j) laserstop(j)],...
%             [min(y(:))-scalefactor*yrange max(y(:))+scalefactor*yrange max(y(:))+scalefactor*yrange min(y(:))-scalefactor*yrange],...
%             pcolor, 'FaceAlpha', 1, 'EdgeColor', 'none');
        
        patch([laserstart(j) laserstart(j) laserstop(j) laserstop(j)],...
            [max(y(:))+scalefactor*yrange max(y(:))+2*scalefactor*yrange max(y(:))+2*scalefactor*yrange max(y(:))+scalefactor*yrange],...
            pcolor, 'FaceAlpha', 1, 'EdgeColor', 'none');
    end
end
hp = zeros(1, size(y, 2));
for j = 1:size(y, 2)
    hp(j) = plot(t, y(:, j), '-', 'Color', lcolor{j});
end
xlim(time_range);
% ylim([min(y(:)) max(y(:))]);
ylim([min(y(:)) max(y(:))+2*scalefactor*yrange]);
xlabel('Time (s)');
ylabel(ylabel_text);
box off;
legend(hp, llabel, 'Location', 'northeast', 'FontSize', 12);
legend('boxoff');
set(gca, 'FontSize', 12, 'ButtonDownFcn', @extract_figure);
end