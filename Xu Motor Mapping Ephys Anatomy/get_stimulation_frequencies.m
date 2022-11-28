%% get stimulation frequencies from stimulation parameter files
curpwd = pwd;
try
   cd(pathname); 
end
[filename, pathname] = uigetfile('*.mat', 'Pick all of the stimulation parameter files', 'MultiSelect', 'on');
if isequal(filename, 0)
    cd(curpwd);
    return;
end

N = numel(filename);
frequency = nan(1, N);

for i = 1:N
    temp = load([pathname filename{i}]);
    frequency(i) = temp.StimulusParameter.Frequency;
end

cd(curpwd);