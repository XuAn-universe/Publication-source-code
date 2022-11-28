%%
FrameRate = 120;
trange = FrameRate*0.5;
stp = 3;
edges = 0:stp:90;
x = stp/2:stp:90-stp/2;

orientation_all_NI = result.orientationxylight_all_NI{1};
orientation_all_I = result.orientationxylight_all_I{1};
orientationbite_all_NI = result.orientationxybite_NI{1};
orientationbite_all_I = result.orientationxybite_I{1};
orientationbite_all_NI = orientationbite_all_NI(trange+1, :);
orientationbite_all_I = orientationbite_all_I(trange+1, :);
[p_NI, ~] = histcounts(orientation_all_NI, edges, 'Normalization', 'probability');
[p_I, ~] = histcounts(orientation_all_I, edges, 'Normalization', 'probability');
[pbite_NI, ~] = histcounts(orientationbite_all_NI, edges, 'Normalization', 'probability');
[pbite_I, ~] = histcounts(orientationbite_all_I, edges, 'Normalization', 'probability');

figure;
plot(x, p_NI, '-k');
hold on;
plot(x, p_I, '-', 'Color', [0 0 0.5]);
plot(x, pbite_NI, '-', 'Color', [0.8 0.8 0.8]);
plot(x, pbite_I, '-b');
xlabel(['Orientation (' char(176) ')']);
ylabel('Probability');
set(gca, 'FontSize', 12, 'TickDir', 'out');
box off;