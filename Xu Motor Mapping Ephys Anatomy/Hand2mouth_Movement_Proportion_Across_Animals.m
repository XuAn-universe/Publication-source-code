%%
a(a(:, 1) == 0, :) = [];
try
    clear result;
end
result(1, :) = [sum(a(:, 2)) sum(a(:, 3) == 1) sum(a(:, 3) == 2)]/size(a, 1);

%%
a(a(:, 1) == 0, :) = [];
result(end+1, :) = [sum(a(:, 2)) sum(a(:, 3) == 1) sum(a(:, 3) == 2)]/size(a, 1);

%% for plotting
figure;
hold on;
for i = 1:size(result, 2)
    bar(i, mean(result(:, i)), 0.6, 'FaceColor', [0.5 0.5 0.5]);
    errorbar(i, mean(result(:, i)), std(result(:, i))/sqrt(numel(result(:, i))), 'k', 'LineStyle', 'none', 'CapSize', 15);
end
for i = 1:size(result, 1)
    plot(result(i, :), '-', 'Color', [0.8 0.8 0.8], 'LineWidth', 1);
end
if size(result, 2) == 3
    ylabel('Hand-to-mouth movement probability');
    set(gca, 'XTick', 1:3, 'XTickLabel', {'C+B', 'C', 'B'}, 'TickDir', 'out', 'FontSize', 12);
elseif size(result, 2) == 2
    ylabel('Head-to-mouth movement probability');
    set(gca, 'XTick', 1:2, 'XTickLabel', {'Pre', 'During'}, 'TickDir', 'out', 'FontSize', 12);
end
xlim([0.5 size(result, 2)+0.5]);

%% for hand-to-mouth movement
curpwd = pwd;
try
   cd(pathname); 
end
[filename, pathname] = uigetfile('*.xlsx', 'Pick all of the data set', 'MultiSelect', 'on');
if isequal(filename, 0)
    cd(curpwd);
    return;
end

N = numel(filename);

ppre = nan(1, N);
pcontra = nan(1, N);
pbi = nan(1, N);

for i = 1:N
    [num, ~, ~] = xlsread([pathname filename{i}]);
    includeID = num(:, 2);
    data = num(includeID == 1, :);
    ppre(i) = sum(data(:, 3) ~= 0)/size(data, 1);
    pcontra(i) = sum(data(:, 4) == 1)/size(data, 1);
    pbi(i) = sum(data(:, 4) == 2)/size(data, 1);
end
result = [ppre' pcontra' pbi'];
cd(curpwd);

%% for head-to-hand movement
curpwd = pwd;
try
   cd(pathname); 
end
[filename, pathname] = uigetfile('*.xlsx', 'Pick all of the data set', 'MultiSelect', 'on');
if isequal(filename, 0)
    cd(curpwd);
    return;
end

N = numel(filename);

ppre = nan(1, N);
pduring = nan(1, N);

for i = 1:N
    [num, ~, ~] = xlsread([pathname filename{i}]);
    includeID = num(:, 2);
    data = num(includeID == 1, :);
    ppre(i) = sum(data(:, 6) ~= 0)/size(data, 1);
    pduring(i) = sum(data(:, 7) == 1)/size(data, 1);
end
result = [ppre' pduring'];
cd(curpwd);