function GCaMP2sit(app, Exp_Path)
Exp_Paths{1} = Exp_Path;

value = app.ListBox.Value;
if ~isempty(value)
    for i = 1:numel(value)
        Exp_Paths{end+1} = app.ListBox.Items{value(i)};
    end
end

pre = app.preEditField.Value;
post = app.postEditField.Value;

zscore_nose_all = [];
zscore_nose_random_all = [];
zscore_pawl_all = [];
zscore_pawl_random_all = [];

zscore_firstbite_all = [];
zscore_firstbite_random_all = [];
firstbite_all = [];

for i = 1:numel(Exp_Paths)
    load([Exp_Paths{i} '\Analysis_Session.mat'], 'Video_annotation');
    
    % get bite events
    try
        audiolocation = Exp_Paths{i}(1:end-7);
        temp = load([audiolocation '\Detected_Bite_Events.mat']);
        Bite_events = temp.Audio_analysis;
    catch
        Bite_events = [];
    end
    
    % get photometry data
    try
        fpdata_all = load([audiolocation '\FPData.mat']);
        nchannel = size(fpdata_all.zsignal_all, 1);
        for j = 1:nchannel
            lgdtext{j} = ['Channel ' num2str(j)];
        end
    catch
        errordlg('Fiber photometry data is missing!', 'Error');
        return;
    end

    for j = 1:numel(Video_annotation)
        if ~Video_annotation(j).Disgard
            if ~isempty(Bite_events)
                bite_timestamps = Bite_events(j).time_bites;
            end
            
            fpdata = fpdata_all.zsignal_all(:, j);
            fpdata_zsignal = [];
            fpdata_t = [];
            for k = 1:nchannel
                fpdata_zsignal(:, k) = fpdata{k}(:, 2);
            end
            fpdata_t = fpdata{1}(:, 1);
            
            if ~isempty(Video_annotation(j).tsit_nose) && Video_annotation(j).tsit_nose ~= 0
                [t, zscore, zscore_random] = AlignSignal2Event(fpdata_t, fpdata_zsignal, Video_annotation(j).tsit_nose, pre, post);
                zscore_nose_all = [zscore_nose_all zscore];
                zscore_nose_random_all = [zscore_nose_random_all zscore_random];
            end
            
            if ~isempty(Video_annotation(j).tsit_pawl) && Video_annotation(j).tsit_pawl ~= 0
                [t, zscore, zscore_random] = AlignSignal2Event(fpdata_t, fpdata_zsignal, Video_annotation(j).tsit_pawl, pre, post);
                zscore_pawl_all = [zscore_pawl_all zscore];
                zscore_pawl_random_all = [zscore_pawl_random_all zscore_random];
            end
            
            if ~isempty(Bite_events)
                [t, zscore, zscore_random] = AlignSignal2Event(fpdata_t, fpdata_zsignal, bite_timestamps(1), pre, post);
                zscore_firstbite_all = [zscore_firstbite_all zscore];
                zscore_firstbite_random_all = [zscore_firstbite_random_all zscore_random];
                
                firstbite_all = [firstbite_all; bite_timestamps(1)];
            end
        end
    end
end

for i = 1:nchannel
    figure;
    if ~isempty(Bite_events)
        subplot(1, 3, 1);
        yyaxis left;
        plot_gcamp_align2event(t, zscore_nose_random_all(:, :, i), zscore_nose_all(:, :, i), ['aligned to nose moving up: ' lgdtext{i}]);
        set(gca, 'YColor', [0 0.5 0]);
        
        yyaxis right;
        h1 = cdfplot(firstbite_all);
        line([min(firstbite_all) min(firstbite_all)], [0 1], 'Color', [1 0 0], 'LineStyle', '--', 'LineWidth', 1);
        set(h1, 'Color', [1 0 0], 'LineWidth', 1.5);
        set(gca, 'GridLineStyle', ':', 'YColor', [1 0 0]);
        xlabel('Time (s)');
        ylabel('Cumulative Probability');
        title(['aligned to nose moving up: ' lgdtext{i}]);
        
        subplot(1, 3, 2);
        yyaxis left;
        plot_gcamp_align2event(t, zscore_pawl_random_all(:, :, i), zscore_pawl_all(:, :, i), ['aligned to left paw moving to food: ' lgdtext{i}]);
        set(gca, 'YColor', [0 0.5 0]);
        
        yyaxis right;
        h1 = cdfplot(firstbite_all);
        line([min(firstbite_all) min(firstbite_all)], [0 1], 'Color', [1 0 0], 'LineStyle', '--', 'LineWidth', 1);
        set(h1, 'Color', [1 0 0], 'LineWidth', 1.5);
        set(gca, 'GridLineStyle', ':', 'YColor', [1 0 0]);
        xlabel('Time (s)');
        ylabel('Cumulative Probability');
        title(['aligned to left paw moving to food: ' lgdtext{i}]);
        
        subplot(1, 3, 3);
        plot_gcamp_align2event(t, zscore_firstbite_random_all(:, :, i), zscore_firstbite_all(:, :, i), ['aligned to first bite: ' lgdtext{i}]);
        set(gca, 'ButtonDownFcn', @extract_figure);
    else
        subplot(1, 2, 1);
        plot_gcamp_align2event(t, zscore_nose_random_all(:, :, i), zscore_nose_all(:, :, i), ['aligned to nose moving up: ' lgdtext{i}]);
        set(gca, 'ButtonDownFcn', @extract_figure);
        
        subplot(1, 2, 2);
        plot_gcamp_align2event(t, zscore_pawl_random_all(:, :, i), zscore_pawl_all(:, :, i), ['aligned to left paw moving to food: ' lgdtext{i}]);
        set(gca, 'ButtonDownFcn', @extract_figure);
    end
end