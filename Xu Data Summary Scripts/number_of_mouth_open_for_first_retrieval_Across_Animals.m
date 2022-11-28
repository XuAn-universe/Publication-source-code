%*---------------------------------------------------------------------*
% Huang Lab
% Duke University
% Author : Xu An, Oct 2021
% xu.an@duke.edu
% Version: 1.0
%*---------------------------------------------------------------------*
%% number of jaw movement for first retrieval
curpwd = pwd;
try
   cd(pathname); 
end
[filename_saline, pathname] = uigetfile('*.mat', 'Pick all of the saline data set', 'MultiSelect', 'on');
if isequal(filename_saline, 0)
    cd(curpwd);
    return;
end
N = numel(filename_saline);
njaw4retrieval1_avg_saline = nan(1, N);
njaw4retrieval1_all_saline = cell(1, N);
njaw4retrieval1_saline = [];
for i = 1:N
    temp = load([pathname filename_saline{i}]);
    result = temp.result;
    temp = result.njaw4retrieval1;
    temp = rmmissing(temp);
    njaw4retrieval1_avg_saline(i) = mean(temp);
    njaw4retrieval1_all_saline{i} = temp;
    njaw4retrieval1_saline = [njaw4retrieval1_saline temp];
end

cd(pathname); 
[filename_muscimol, pathname] = uigetfile('*.mat', 'Pick all of the muscimol data set', 'MultiSelect', 'on');
if isequal(filename_muscimol, 0)
    return;
end
N = numel(filename_muscimol);
njaw4retrieval1_avg_muscimol = nan(1, N);
njaw4retrieval1_all_muscimol = cell(1, N);
njaw4retrieval1_muscimol = [];
for i = 1:N
    temp = load([pathname filename_muscimol{i}]);
    result = temp.result;
    temp = result.njaw4retrieval1;
    temp = rmmissing(temp);
    njaw4retrieval1_avg_muscimol(i) = mean(temp);
    njaw4retrieval1_all_muscimol{i} = temp;
    njaw4retrieval1_muscimol = [njaw4retrieval1_muscimol temp];
end

resampling = 10000;
figure;
permutation_test(njaw4retrieval1_all_saline, njaw4retrieval1_all_muscimol, 'both', resampling, 'difference of # of mouth open');
figure;
bootstrapping_test(njaw4retrieval1_all_saline, njaw4retrieval1_all_muscimol, 'both', resampling, 'difference of # of mouth open');

figure;
paired_plot(njaw4retrieval1_avg_saline, njaw4retrieval1_avg_muscimol, '# of mouth open', filename_saline);

compare2distributions(njaw4retrieval1_saline, njaw4retrieval1_muscimol, '# of mouth open', '');

violin_boxplot(njaw4retrieval1_saline, njaw4retrieval1_muscimol, '# of mouth open');

cd(curpwd);

%% start time of successful retrieval
curpwd = pwd;
try
   cd(pathname); 
end
[filename_saline, pathname] = uigetfile('*.mat', 'Pick all of the saline data set', 'MultiSelect', 'on');
if isequal(filename_saline, 0)
    cd(curpwd);
    return;
end
N = numel(filename_saline);
RetrievalStart1st_avg_saline = nan(1, N);
RetrievalStart1st_all_saline = cell(1, N);
RetrievalStart1st_saline = [];
for i = 1:N
    temp = load([pathname filename_saline{i}]);
    result = temp.result;
    temp = result.RetrievalStart1st;
    temp = rmmissing(temp);
    RetrievalStart1st_avg_saline(i) = mean(temp);
    RetrievalStart1st_all_saline{i} = temp;
    RetrievalStart1st_saline = [RetrievalStart1st_saline temp];
end

cd(pathname); 
[filename_muscimol, pathname] = uigetfile('*.mat', 'Pick all of the muscimol data set', 'MultiSelect', 'on');
if isequal(filename_muscimol, 0)
    return;
end
N = numel(filename_muscimol);
RetrievalStart1st_avg_muscimol = nan(1, N);
RetrievalStart1st_all_muscimol = cell(1, N);
RetrievalStart1st_muscimol = [];
for i = 1:N
    temp = load([pathname filename_muscimol{i}]);
    result = temp.result;
    temp = result.RetrievalStart1st;
    temp = rmmissing(temp);
    RetrievalStart1st_avg_muscimol(i) = mean(temp);
    RetrievalStart1st_all_muscimol{i} = temp;
    RetrievalStart1st_muscimol = [RetrievalStart1st_muscimol temp];
end

resampling = 10000;
figure;
permutation_test(RetrievalStart1st_all_saline, RetrievalStart1st_all_muscimol, 'both', resampling, 'difference of successful retrieval start');
figure;
bootstrapping_test(RetrievalStart1st_all_saline, RetrievalStart1st_all_muscimol, 'both', resampling, 'difference of successful retrieval start');

figure;
paired_plot(RetrievalStart1st_avg_saline, RetrievalStart1st_avg_muscimol, 'Successful retrieval start (s)', filename_saline);

compare2distributions(RetrievalStart1st_saline, RetrievalStart1st_muscimol, 'Successful retrieval start (s)', '');

violin_boxplot(RetrievalStart1st_saline, RetrievalStart1st_muscimol, 'Successful retrieval start (s)');

cd(curpwd);

%% pasta detection time
optogenetic = 0; % condition between muscimol and optogenetics
if ~optogenetic
    curpwd = pwd;
    try
        cd(pathname);
    end
    [filename_saline, pathname] = uigetfile('*.mat', 'Pick all of the saline data set', 'MultiSelect', 'on');
    if isequal(filename_saline, 0)
        cd(curpwd);
        return;
    end
    N = numel(filename_saline);
    timedt_avg_saline = nan(1, N);
    timedt_all_saline = cell(1, N);
    timedt_saline = [];
    for i = 1:N
        temp = load([pathname filename_saline{i}]);
        timedt = temp.timedt;
        timedt_avg_saline(i) = mean(timedt(:, 1));
        timedt_all_saline{i} = timedt(:, 1)';
        timedt_saline = [timedt_saline timedt(:, 1)'];
    end
    
    cd(pathname);
    [filename_muscimol, pathname] = uigetfile('*.mat', 'Pick all of the muscimol data set', 'MultiSelect', 'on');
    if isequal(filename_muscimol, 0)
        return;
    end
    N = numel(filename_muscimol);
    timedt_avg_muscimol = nan(1, N);
    timedt_all_muscimol = cell(1, N);
    timedt_muscimol = [];
    for i = 1:N
        temp = load([pathname filename_muscimol{i}]);
        timedt = temp.timedt;
        timedt_avg_muscimol(i) = mean(timedt(:, 1));
        timedt_all_muscimol{i} = timedt(:, 1)';
        timedt_muscimol = [timedt_muscimol timedt(:, 1)'];
    end
else
    curpwd = pwd;
    try
        cd(pathname);
    end
    [filename, pathname] = uigetfile('*.mat', 'Pick all of the data set', 'MultiSelect', 'on');
    if isequal(filename, 0)
        return;
    end
    N = numel(filename);
    timedt_avg_NI = nan(1, N);
    timedt_all_NI = cell(1, N);
    timedt_NI = [];
    timedt_avg_I = nan(1, N);
    timedt_all_I = cell(1, N);
    timedt_I = [];
    for i = 1:N
        temp = load([pathname filename{i}]);
        timedt = temp.timedt;
        timedt_avg_NI(i) = mean(timedt(timedt(:, 2)==0, 1));
        timedt_all_NI{i} = timedt(timedt(:, 2)==0, 1)';
        timedt_NI = [timedt_NI timedt(timedt(:, 2)==0, 1)'];
        timedt_avg_I(i) = mean(timedt(timedt(:, 2)==1, 1));
        timedt_all_I{i} = timedt(timedt(:, 2)==1, 1)';
        timedt_I = [timedt_I timedt(timedt(:, 2)==1, 1)'];
    end
    timedt_all_saline = timedt_all_NI;
    timedt_all_muscimol = timedt_all_I;
    timedt_avg_saline = timedt_avg_NI;
    timedt_avg_muscimol = timedt_avg_I;
    timedt_saline = timedt_NI;
    timedt_muscimol = timedt_I;
    filename_saline = filename;
end

resampling = 10000;
figure;
permutation_test(timedt_all_saline, timedt_all_muscimol, 'both', resampling, 'difference of pasta detection time');
figure;
bootstrapping_test(timedt_all_saline, timedt_all_muscimol, 'both', resampling, 'difference of pasta detection time');

figure;
paired_plot(timedt_avg_saline, timedt_avg_muscimol, 'Pasta detection (s)', filename_saline);

compare2distributions(timedt_saline, timedt_muscimol, 'Pasta detection (s)', '');

violin_boxplot(timedt_saline, timedt_muscimol, 'Pasta detection (s)');

cd(curpwd);