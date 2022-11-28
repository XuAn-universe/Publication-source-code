%*---------------------------------------------------------------------*
% Huang Lab
% Duke University
% Author : Xu An, Nov 2021
% xu.an@duke.edu
% Version: 1.0
%*---------------------------------------------------------------------*
%% pasta orientation change analysis across mice
curpwd = pwd;
try
   cd(pathname); 
end
[filename, pathname] = uigetfile('*.mat', 'Pick all of the data set', 'MultiSelect', 'on');
if isequal(filename, 0)
    cd(curpwd);
    return;
end

FrameRate = 120;
trange = 0.5;
trange = trange*FrameRate;
t = (-trange:1:trange)'/FrameRate;

N = numel(filename);
orientationchange_NI = nan(1, N);
orientationchange_I = nan(1, N);
orientationchangebite_NI = nan(2*trange+1, N);
orientationchangebite_I = nan(2*trange+1, N);
for i = 1:N
    temp = load([pathname filename{i}]);
    result = temp.result;
    temp = result.orientationxychangelight_mean_NI;
    orientationchange_NI(i) = temp;
    temp = result.orientationxychangelight_mean_I;
    orientationchange_I(i) = temp;
    
    orientationchangebite_NI(:, i) = mean(result.orientationxychangebite_NI{1}, 2, 'omitnan');
    if ~isempty(result.orientationxychangebite_I{1})
        orientationchangebite_I(:, i) = mean(result.orientationxychangebite_I{1}, 2, 'omitnan');
    end
end
figure;
paired_plot(orientationchange_NI, orientationchange_I, ['Orientation (' char(176) '/s)'], filename, 'bar');
cd(curpwd);

figure;
plot_tj_individuals(t, orientationchangebite_NI, [0.8 0.8 0.8], [0 0 0], 'Time (s)', ['Orientation (' char(176) '/s)'], '');
plot_tj_individuals(t, orientationchangebite_I, [0 1 1], [0 0 1], 'Time (s)', ['Orientation (' char(176) '/s)'], '');
figure;
plot_tj_acrosstrial(t, orientationchangebite_NI, [0 0 0], [0 0 0], ['Orientation (' char(176) '/s)']);
plot_tj_acrosstrial(t, orientationchangebite_I, [0 0 1], [0 0 1], ['Orientation (' char(176) '/s)']);