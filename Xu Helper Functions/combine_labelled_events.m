%%
n = 2;
switch n
    case 2
        trials_NI = [trials_NI1 trials_NI2];
        trials_I = [trials_I1 trials_I2];
    case 3
        trials_NI = [trials_NI1 trials_NI2 trials_NI3];
        trials_I = [trials_I1 trials_I2 trials_I3];
    case 4
        trials_NI = [trials_NI1 trials_NI2 trials_NI3 trials_NI4];
        trials_I = [trials_I1 trials_I2 trials_I3 trials_I4];    
end

curpwd = pwd;
try
    cd(path);
end
[file, path] = uiputfile('*.mat', 'Save result');
if file ~= 0
    save([path file], 'trials_NI', 'trials_I');
end
cd(curpwd);