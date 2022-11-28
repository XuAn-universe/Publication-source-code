function AverageOmission(app, Exp_Path)
Exp_Paths{1} = Exp_Path;

value = app.ListBox.Value;
if ~isempty(value)
    for i = 1:numel(value)
        Exp_Paths{end+1} = app.ListBox.Items{value(i)};
    end
end

pre = 4;
post = app.postEditField.Value;

zscore_all = cell(0);
SampleRate = [];
nomission = 0;

for i = 1:numel(Exp_Paths)
    load([Exp_Paths{i} '\Analysis_Session.mat'], 'Video_annotation');
    
    % get photometry data
    try
        signallocation = Exp_Paths{i}(1:end-7);
        fpdata_all = load([signallocation '\FPData.mat']);
        nchannel = size(fpdata_all.zsignal_all, 1);
        for j = 1:nchannel
            lgdtext{j} = ['Channel ' num2str(j)];
        end
    catch
        errordlg('Fiber photometry data is missing!', 'Error');
        return;
    end

    for j = 1:numel(Video_annotation)
        if ~Video_annotation(j).Disgard && Video_annotation(j).Omission
            nomission = nomission+1;
            fpdata = fpdata_all.zsignal_all(:, j);
            fpdata_zsignal = [];
            fpdata_t = [];
            for k = 1:nchannel
                fpdata_zsignal(:, k) = fpdata{k}(:, 2);
            end
            fpdata_t = fpdata{1}(:, 1);
            
            zscore_all{nomission} = fpdata_zsignal;
            SampleRate(nomission) = 1/mean(diff(fpdata_t));
        end
    end
end

SampleRate = mean(SampleRate);
t = -pre:1/SampleRate:post;
for i = 1:nchannel
    zscore_matrix = nan(numel(t), nomission);
    for j = 1:nomission
        ntemp = size(zscore_all{j}, 1);
        if ntemp >= numel(t)
            zscore_matrix(:, j) = zscore_all{j}(1:numel(t), i);
        else
            zscore_matrix(1:ntemp, j) = zscore_all{j}(:, i);
        end
    end
    figure;
    plot_tj_individuals(t, zscore_matrix, [0.75 0.75 0.75], [0 0 0], 'Time (s)', 'Z score', lgdtext{i});
end