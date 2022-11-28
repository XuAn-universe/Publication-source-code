function PCAZaxis(app, Exp_Path, FrameRate, nmedian)
Exp_Paths{1} = Exp_Path;

value = app.ListBox.Value;
if ~isempty(value)
    for i = 1:numel(value)
        Exp_Paths{end+1} = app.ListBox.Items{value(i)};
    end
end

if app.TrackingDataCheckBox.Value
    label_table = table2array(app.UITable.Data);
else
    helpdlg('''Tracking Data'' needs to be checked');
    return;
end
pointIDs = find(label_table(:, 2));
nPC = numel(pointIDs);
nPC_plot = app.PCEditField.Value;

if app.sOnTwiceButton.Value
    ninhibition = 2;
    tinhibition = [4 13; 8 17];
elseif app.NoLightButton.Value || app.WholeTrialOnButton.Value
    ninhibition = 1;
    tinhibition = [0 inf];
elseif app.sDelay4sOnButton.Value
    ninhibition = 1;
    tinhibition = [4 8];
elseif app.sDelay4sOnButton_2.Value
    ninhibition = 1;
    tinhibition = [0 4];
end

laserstart_all = [];
laserstop_all = [];
Zmatrix = [];
n_all = 0;
FrameRate = round(FrameRate);

bite_timestamps_all = cell(0);
laser_timestamps_all = cell(0);

for i = 1:numel(Exp_Paths)
    load([Exp_Paths{i} '\Analysis_Session.mat'], 'Video_annotation');
    try
        audiolocation = Exp_Paths{i}(1:end-7);
        temp = load([audiolocation '\Detected_Bite_Events.mat']);
        Bite_events = temp.Audio_analysis;
    catch
        errordlg('Please detect the bites first.', 'Error');
        return;
    end
    for j = 1:numel(Video_annotation)
        if ~Video_annotation(j).Disgard
            n_all = n_all+1;
            Zsingletrial = [];
            for ii = 1:nPC
                [~, ~, ~, ~, ~, ~, ~, z, ~, ~, laserstart, laserstop] = trajectory_postprocessing(pointIDs(ii), Exp_Paths{i}, j,...
                    label_table, nmedian, FrameRate);
                Zsingletrial = [Zsingletrial z(:, 1)];
            end
            t = (1:size(Zsingletrial, 1))'/FrameRate;
            
            bite_timestamps = Bite_events(j).time_bites;
            bite_timestamps_all{n_all} = bite_timestamps;
            if app.NoLightButton.Value
                laserstart = [];
                laserstop = [];
            else
                laser_timestamps = Bite_events(j).laser_timestamps;
                laser_timestamps_all{n_all} = laser_timestamps;
            end
            
            Zsingletrial(t < bite_timestamps(1) | t > bite_timestamps(end), :) = NaN;
            framecounter(n_all) = size(Zsingletrial, 1);
            if ~isempty(laserstart)
                trialtype(n_all) = 1;
                if numel(laserstart) ~= ninhibition
                    laserstart_all = [laserstart_all round((Bite_events(j).laser_timestamps(1, :))'*FrameRate)];
                    laserstop_all = [laserstop_all round((Bite_events(j).laser_timestamps(2, :))'*FrameRate)];
                else
                    laserstart_all = [laserstart_all laserstart];
                    laserstop_all = [laserstop_all laserstop];
                end
            else
                trialtype(n_all) = 0;
            end
            Zmatrix = [Zmatrix; Zsingletrial];
        end
    end
end

% a = rand(1000, 3);
% a = bsxfun(@minus, a, mean(a));
% [coeff, score, latent, tsquared, explained] = pca(a);
% a1 = score(:, 1)*coeff(:, 1)';
% a2 = score(:, 2)*coeff(:, 2)';
% a3 = score(:, 3)*coeff(:, 3)';
% b = a1+a2+a3;
% a12 = score(:, 1:2)*coeff(:, 1:2)';
% c = a12+a3;
% a, b, c are equal
[coeff, score, latent, tsquared, explained] = pca(Zmatrix);
assignin('base', 'variance_explained', explained);
assignin('base', 'coeff', coeff);

n_I = sum(trialtype);
n_NI = n_all-n_I;
frames_I = max(framecounter(trialtype == 1));
frames_NI = max(framecounter(trialtype == 0));

score_all_NI = cell(0);
speed_all_NI = cell(0);
acceleration_all_NI = cell(0);

score_all_I = cell(0);
speed_all_I = cell(0);
acceleration_all_I = cell(0);

scorebite_NI = [];
speedbite_NI = [];
accelerationbite_NI = [];

scorebite_I = [];
speedbite_I = [];
accelerationbite_I = [];

nbite_NI = zeros(1, ninhibition);
nbite_I = zeros(1, ninhibition);

for i = 1:n_all
    if i == 1
        score_singletrial = score(1:framecounter(1), :);
    else
        score_singletrial = score(sum(framecounter(1:i-1))+1:sum(framecounter(1:i)), :);
    end
    speed_singletrial = diff(score_singletrial);
    speed_singletrial = [NaN(1, nPC); speed_singletrial];
    acceleration_singletrial = diff(speed_singletrial);
    acceleration_singletrial = [NaN(1, nPC); acceleration_singletrial];
    
    t = (1:framecounter(i))'/FrameRate;
    
    if trialtype(i)
        laserstart = laserstart_all(:, sum(trialtype(1:i)));
        laserstop = laserstop_all(:, sum(trialtype(1:i)));
    else
        laserstart = [];
        laserstop = [];
    end
    
    bite_timestamps = bite_timestamps_all{i};
    if ~app.NoLightButton.Value
        laser_timestamps = laser_timestamps_all{i};
    end
    
    if ~isempty(laserstart)
        score_all_I{end+1} = score_singletrial;
        speed_all_I{end+1} = speed_singletrial;
        acceleration_all_I{end+1} = acceleration_singletrial;
        
        if ~isempty(bite_timestamps)
            for k = 1:ninhibition
                for m = 1:numel(bite_timestamps)
                    timestamp = bite_timestamps(m);
                    if timestamp >= laser_timestamps(2*k-1) && timestamp <= laser_timestamps(2*k)
                        [~, id] = min(abs(t-timestamp));
                        if id-FrameRate >= 1
                            nbite_I(k) = nbite_I(k)+1;
                            scorebite_I(:, :, nbite_I(k), k) = score_singletrial(id-FrameRate:id+FrameRate, :);
                            speedbite_I(:, :, nbite_I(k), k) = abs(speed_singletrial(id-FrameRate:id+FrameRate, :));
                            accelerationbite_I(:, :, nbite_I(k), k) = abs(acceleration_singletrial(id-FrameRate:id+FrameRate, :));
                        end
                    end
                end
            end
        end
    else
        score_all_NI{end+1} = score_singletrial;
        speed_all_NI{end+1} = speed_singletrial;
        acceleration_all_NI{end+1} = acceleration_singletrial;
        
        if ~isempty(bite_timestamps)
            for k = 1:ninhibition
                for m = 1:numel(bite_timestamps)
                    timestamp = bite_timestamps(m);
                    if timestamp >= tinhibition(2*k-1) && timestamp <= tinhibition(2*k)
                        [~, id] = min(abs(t-timestamp));
                        if id-FrameRate >= 1
                            nbite_NI(k) = nbite_NI(k)+1;
                            scorebite_NI(:, :, nbite_NI(k), k) = score_singletrial(id-FrameRate:id+FrameRate, :);
                            speedbite_NI(:, :, nbite_NI(k), k) = abs(speed_singletrial(id-FrameRate:id+FrameRate, :));
                            accelerationbite_NI(:, :, nbite_NI(k), k) = abs(acceleration_singletrial(id-FrameRate:id+FrameRate, :));
                        end
                    end
                end
            end
        end
    end
end

figure;
hold on;
for i = 1:n_NI
    plot(score_all_NI{i}(:, 1), score_all_NI{i}(:, 2), 'xk');
end
figure;
hold on;
for i = 1:n_I
    plot(score_all_I{i}(:, 1), score_all_I{i}(:, 2), 'xk');
end

score_matrix_NI = NaN(frames_NI, nPC, n_NI);
speed_matrix_NI = NaN(frames_NI, nPC, n_NI);
acceleration_matrix_NI = NaN(frames_NI, nPC, n_NI);

for i = 1:ninhibition
    scorelight_all_NI{i} = [];
    speedlight_all_NI{i} = [];
    accelerationlight_all_NI{i} = [];
    
    scorelight_mean_NI{i} = [];
    speedlight_mean_NI{i} = [];
    accelerationlight_mean_NI{i} = [];
    
    scorelight_var_NI{i} = [];
    speedlight_var_NI{i} = [];
    accelerationlight_var_NI{i} = [];
end
for i = 1:n_NI
    frames_temp = size(score_all_NI{i}, 1);
    score_matrix_NI(1:frames_temp, :, i) = score_all_NI{i};
    speed_matrix_NI(1:frames_temp, :, i) = abs(speed_all_NI{i});
    acceleration_matrix_NI(1:frames_temp, :, i) = abs(acceleration_all_NI{i});
    
    t = (1:frames_temp)'/FrameRate;
    for j = 1:ninhibition
        scorelight_all_NI{j} = [scorelight_all_NI{j}; score_all_NI{i}(t >= tinhibition(2*j-1) & t <= tinhibition(2*j), :)];
        speedlight_all_NI{j} = [speedlight_all_NI{j}; abs(speed_all_NI{i}(t >= tinhibition(2*j-1) & t <= tinhibition(2*j), :))];
        accelerationlight_all_NI{j} = [accelerationlight_all_NI{j}; abs(acceleration_all_NI{i}(t >= tinhibition(2*j-1) & t <= tinhibition(2*j), :))];
        
        scorelight_mean_NI{j} = [scorelight_mean_NI{j}; mean(score_all_NI{i}(t >= tinhibition(2*j-1) & t <= tinhibition(2*j), :), 'omitnan')];
        speedlight_mean_NI{j} = [speedlight_mean_NI{j}; mean(abs(speed_all_NI{i}(t >= tinhibition(2*j-1) & t <= tinhibition(2*j), :)), 'omitnan')];
        accelerationlight_mean_NI{j} = [accelerationlight_mean_NI{j}; mean(abs(acceleration_all_NI{i}(t >= tinhibition(2*j-1) & t <= tinhibition(2*j), :)), 'omitnan')];
        
        scorelight_var_NI{j} = [scorelight_var_NI{j}; var(score_all_NI{i}(t >= tinhibition(2*j-1) & t <= tinhibition(2*j), :), 'omitnan')];
        speedlight_var_NI{j} = [speedlight_var_NI{j}; var(abs(speed_all_NI{i}(t >= tinhibition(2*j-1) & t <= tinhibition(2*j), :)), 'omitnan')];
        accelerationlight_var_NI{j} = [accelerationlight_var_NI{j}; var(abs(acceleration_all_NI{i}(t >= tinhibition(2*j-1) & t <= tinhibition(2*j), :)), 'omitnan')];
        
        if i == 1
            for k = 1:nPC_plot
                figure;
                haxes(j, k) = axes;
                hold on;
                xlabel('Time (s)');
                ylabel('Principal Posistion (mm)');
                title(['PM ' num2str(k)]);
                set(gca, 'XLim', [tinhibition(2*j-1) tinhibition(2*j)], 'TickLength', [0 0], 'FontSize', 12);
                box off;
            end
            figure;
            haxes(j, k+1) = axes;
            hold on;
            xlabel('PM 1');
            ylabel('PM 2');
            title(['Inhibition ' num2str(j)]);
            set(gca, 'TickLength', [0 0], 'FontSize', 12);
            box off;
            figure;
            haxes(j, k+2) = axes;
            hold on;
            xlabel('PM 1');
            ylabel('PM 3');
            title(['Inhibition ' num2str(j)]);
            set(gca, 'TickLength', [0 0], 'FontSize', 12);
            box off;
            figure;
            haxes(j, k+3) = axes;
            hold on;
            xlabel('PM 2');
            ylabel('PM 3');
            title(['Inhibition ' num2str(j)]);
            set(gca, 'TickLength', [0 0], 'FontSize', 12);
            box off;
        end
        for k = 1:nPC_plot
            plot(haxes(j, k), t(t >= tinhibition(2*j-1) & t <= tinhibition(2*j)), score_all_NI{i}(t >= tinhibition(2*j-1) & t <= tinhibition(2*j), k), '-k');
        end
        plot(haxes(j, k+1), score_all_NI{i}(t >= tinhibition(2*j-1) & t <= tinhibition(2*j), 1), score_all_NI{i}(t >= tinhibition(2*j-1) & t <= tinhibition(2*j), 2), 'ok');
        plot(haxes(j, k+2), score_all_NI{i}(t >= tinhibition(2*j-1) & t <= tinhibition(2*j), 1), score_all_NI{i}(t >= tinhibition(2*j-1) & t <= tinhibition(2*j), 3), 'ok');
        plot(haxes(j, k+3), score_all_NI{i}(t >= tinhibition(2*j-1) & t <= tinhibition(2*j), 2), score_all_NI{i}(t >= tinhibition(2*j-1) & t <= tinhibition(2*j), 3), 'ok');
    end
end
assignin('base', 'scorelight_var_NI', scorelight_var_NI);
t_NI = (1:frames_NI)'/FrameRate;
time_range = [app.timerangesEditField.Value app.toEditField.Value];

if ~app.NoLightButton.Value
    score_matrix_I = NaN(frames_I, nPC, n_I);
    speed_matrix_I = NaN(frames_I, nPC, n_I);
    acceleration_matrix_I = NaN(frames_I, nPC, n_I);
    
    for i = 1:ninhibition
        scorelight_all_I{i} = [];
        speedlight_all_I{i} = [];
        accelerationlight_all_I{i} = [];
        
        scorelight_mean_I{i} = [];
        speedlight_mean_I{i} = [];
        accelerationlight_mean_I{i} = [];
        
        scorelight_var_I{i} = [];
        speedlight_var_I{i} = [];
        accelerationlight_var_I{i} = [];
    end
    for i = 1:n_I
        frames_temp = size(score_all_I{i}, 1);
        score_matrix_I(1:frames_temp, :, i) = score_all_I{i};
        speed_matrix_I(1:frames_temp, :, i) = abs(speed_all_I{i});
        acceleration_matrix_I(1:frames_temp, :, i) = abs(acceleration_all_I{i});
        
        for j = 1:ninhibition
            scorelight_all_I{j} = [scorelight_all_I{j}; score_all_I{i}(laserstart_all(j, i):laserstop_all(j, i), :)];
            speedlight_all_I{j} = [speedlight_all_I{j}; abs(speed_all_I{i}(laserstart_all(j, i):laserstop_all(j, i), :))];
            accelerationlight_all_I{j} = [accelerationlight_all_I{j}; abs(acceleration_all_I{i}(laserstart_all(j, i):laserstop_all(j, i), :))];
            
            scorelight_mean_I{j} = [scorelight_mean_I{j}; mean(score_all_I{i}(laserstart_all(j, i):laserstop_all(j, i), :), 'omitnan')];
            speedlight_mean_I{j} = [speedlight_mean_I{j}; mean(abs(speed_all_I{i}(laserstart_all(j, i):laserstop_all(j, i), :)), 'omitnan')];
            accelerationlight_mean_I{j} = [accelerationlight_mean_I{j}; mean(abs(acceleration_all_I{i}(laserstart_all(j, i):laserstop_all(j, i), :)), 'omitnan')];
            
            scorelight_var_I{j} = [scorelight_var_I{j}; var(score_all_I{i}(laserstart_all(j, i):laserstop_all(j, i), :), 'omitnan')];
            speedlight_var_I{j} = [speedlight_var_I{j}; var(abs(speed_all_I{i}(laserstart_all(j, i):laserstop_all(j, i), :)), 'omitnan')];
            accelerationlight_var_I{j} = [accelerationlight_var_I{j}; var(abs(acceleration_all_I{i}(laserstart_all(j, i):laserstop_all(j, i), :)), 'omitnan')];
            
            for k = 1:nPC_plot
                plot(haxes(j, k), (laserstart_all(j, i):laserstop_all(j, i))/FrameRate, score_all_I{i}(laserstart_all(j, i):laserstop_all(j, i), k), '-g');
            end
            plot(haxes(j, k+1), score_all_I{i}(laserstart_all(j, i):laserstop_all(j, i), 1), score_all_I{i}(laserstart_all(j, i):laserstop_all(j, i), 2), 'og');
            plot(haxes(j, k+2), score_all_I{i}(laserstart_all(j, i):laserstop_all(j, i), 1), score_all_I{i}(laserstart_all(j, i):laserstop_all(j, i), 3), 'og');
            plot(haxes(j, k+3), score_all_I{i}(laserstart_all(j, i):laserstop_all(j, i), 2), score_all_I{i}(laserstart_all(j, i):laserstop_all(j, i), 3), 'og');
        end
    end
    laserstart = mean(laserstart_all, 2)/FrameRate;
    laserstop = mean(laserstop_all, 2)/FrameRate;
    t_I = (1:frames_I)'/FrameRate;
    
    for i = 1:nPC_plot
        figure;
        plot_tj_multitrial(time_range, t_NI, squeeze(score_matrix_NI(:, i, :)), t_I, squeeze(score_matrix_I(:, i, :)), laserstart, laserstop, 'Principal Position (mm)', ['PM ' num2str(i)]);
        figure;
        plot_tj_multitrial(time_range, t_NI, squeeze(speed_matrix_NI(:, i, :)), t_I, squeeze(speed_matrix_I(:, i, :)),...
            laserstart, laserstop, 'Speed (mm/s)', ['PM ' num2str(i)]);
        figure;
        plot_tj_multitrial(time_range, t_NI, squeeze(acceleration_matrix_NI(:, i, :)), t_I, squeeze(acceleration_matrix_I(:, i, :)),...
            laserstart, laserstop, 'Absolute Acceleration (mm/s^{2})', ['PM ' num2str(i)]);
    end
else
    for i = 1:nPC_plot
        figure;
        plot_tj_multitrial(time_range, t_NI, score_matrix_NI(:, i, :), [], [], [], [], 'Principal Position (mm)', ['PM ' num2str(i)]);
        figure;
        plot_tj_multitrial(time_range, t_NI, squeeze(speed_matrix_NI(:, i, :)), [], [],...
            [], [], 'Speed (mm/s)', ['PM ' num2str(i)]);
        figure;
        plot_tj_multitrial(time_range, t_NI, squeeze(acceleration_matrix_NI(:, i, :)), [], [],...
            [], [], 'Absolute Acceleration (mm/s^{2})', ['PM ' num2str(i)]);
    end
end
assignin('base', 'scorelight_var_I', scorelight_var_I);

t = (-FrameRate:1:FrameRate)'/FrameRate;
time_range = [t(1) t(end)];
if ~app.NoLightButton.Value
    for i = 1:ninhibition
        for j = 1:nPC_plot
            % trajectories aligned to bite time
            figure;
            plot_tj_multitrial(time_range, t, squeeze(scorebite_NI(:, j, 1:nbite_NI(i), i)), t, squeeze(scorebite_I(:, j, 1:nbite_I(i), i)), 0, NaN, 'Principal Position (mm)', ['PM ' num2str(j) ' (Inhibition #' num2str(i) ')']);
            figure;
            plot_tj_multitrial(time_range, t, squeeze(speedbite_NI(:, j, 1:nbite_NI(i), i)), t, squeeze(speedbite_I(:, j, 1:nbite_I(i), i)), 0, NaN, 'Speed (mm/s)', ['PM ' num2str(j) ' (Inhibition #' num2str(i) ')']);
            figure;
            plot_tj_multitrial(time_range, t, squeeze(accelerationbite_NI(:, j, 1:nbite_NI(i), i)), t, squeeze(accelerationbite_I(:, j, 1:nbite_I(i), i)), 0, NaN, 'Absolute Acceleration (mm/s^{2})', ['PM ' num2str(j) ' (Inhibition #' num2str(i) ')']);
            
            % histograms
            histogram_multitrial(scorelight_all_NI{i}(:, j), scorebite_NI(FrameRate+1, j, 1:nbite_NI(i), i), scorelight_all_I{i}(:, j), scorebite_I(FrameRate+1, j, 1:nbite_I(i), i), 'Principal Position (mm)', ['PM ' num2str(j) ' (Inhibition #' num2str(i) ')']);
            histogram_multitrial(speedlight_all_NI{i}(:, j), speedbite_NI(FrameRate+1, j, 1:nbite_NI(i), i), speedlight_all_I{i}(:, j), speedbite_I(FrameRate+1, j, 1:nbite_I(i), i), 'Speed (mm/s)', ['PM ' num2str(j) ' (Inhibition #' num2str(i) ')']);
            histogram_multitrial(accelerationlight_all_NI{i}(:, j), accelerationbite_NI(FrameRate+1, j, 1:nbite_NI(i), i), accelerationlight_all_I{i}(:, j), accelerationbite_I(FrameRate+1, j, 1:nbite_I(i), i), 'Absolute Acceleration (mm/s^{2})', ['PM ' num2str(j) ' (Inhibition #' num2str(i) ')']);
            
            histogram_multitrial(log(scorelight_var_NI{i}(:, j)), [], log(scorelight_var_I{i}(:, j)), [], 'Log Principal Position Variance (mm^{2})', ['PM ' num2str(j) ' (Inhibition #' num2str(i) ')']);
            histogram_multitrial(log(speedlight_var_NI{i}(:, j)), [], log(speedlight_var_I{i}(:, j)), [], 'Log Speed Variance (mm^{2}/s^{2})', ['PM ' num2str(j) ' (Inhibition #' num2str(i) ')']);
            histogram_multitrial(log(accelerationlight_var_NI{i}(:, j)), [], log(accelerationlight_var_I{i}(:, j)), [], 'Log Absolute Acceleration Variance (mm^{2}/s^{4})', ['PM ' num2str(j) ' (Inhibition #' num2str(i) ')']);
        end
        
        result.scorelight_mean_NI(i, :) = mean(scorelight_mean_NI{i}, 'omitnan');
        result.speedlight_mean_NI(i, :) = mean(speedlight_mean_NI{i}, 'omitnan');
        result.accelerationlight_mean_NI(i, :) = mean(accelerationlight_mean_NI{i}, 'omitnan');
        
        result.scorelight_var_NI(i, :) = mean(scorelight_var_NI{i}, 'omitnan');
        result.speedlight_var_NI(i, :) = mean(speedlight_var_NI{i}, 'omitnan');
        result.accelerationlight_var_NI(i, :) = mean(accelerationlight_var_NI{i}, 'omitnan');
        
        result.scorelight_mean_I(i, :) = mean(scorelight_mean_I{i}, 'omitnan');
        result.speedlight_mean_I(i, :) = mean(speedlight_mean_I{i}, 'omitnan');
        result.accelerationlight_mean_I(i, :) = mean(accelerationlight_mean_I{i}, 'omitnan');
        
        result.scorelight_var_I(i, :) = mean(scorelight_var_I{i}, 'omitnan');
        result.speedlight_var_I(i, :) = mean(speedlight_var_I{i}, 'omitnan');
        result.accelerationlight_var_I(i, :) = mean(accelerationlight_var_I{i}, 'omitnan');
    end
else
    for i = 1:ninhibition
        for j = 1:nPC_plot
            % trajectories aligned to bite time
            figure;
            plot_tj_multitrial(time_range, t, squeeze(scorebite_NI(:, j, 1:nbite_NI(i), i)), [], [], 0, NaN, 'Principal Position (mm)', ['PM ' num2str(j) ' (Inhibition #' num2str(i) ')']);
            figure;
            plot_tj_multitrial(time_range, t, squeeze(speedbite_NI(:, j, 1:nbite_NI(i), i)), [], [], 0, NaN, 'Speed (mm/s)', ['PM ' num2str(j) ' (Inhibition #' num2str(i) ')']);
            figure;
            plot_tj_multitrial(time_range, t, squeeze(accelerationbite_NI(:, j, 1:nbite_NI(i), i)), [], [], 0, NaN, 'Absolute Acceleration (mm/s^{2})', ['PM ' num2str(j) ' (Inhibition #' num2str(i) ')']);
            
            % histograms
            histogram_multitrial(scorelight_all_NI{i}(:, j), scorebite_NI(FrameRate+1, j, 1:nbite_NI(i), i), [], [], 'Principal Position (mm)', ['PM ' num2str(j) ' (Inhibition #' num2str(i) ')']);
            histogram_multitrial(speedlight_all_NI{i}(:, j), speedbite_NI(FrameRate+1, j, 1:nbite_NI(i), i), [], [], 'Speed (mm/s)', ['PM ' num2str(j) ' (Inhibition #' num2str(i) ')']);
            histogram_multitrial(accelerationlight_all_NI{i}(:, j), accelerationbite_NI(FrameRate+1, j, 1:nbite_NI(i), i), [], [], 'Absolute Acceleration (mm/s^{2})', ['PM ' num2str(j) ' (Inhibition #' num2str(i) ')']);
            
            histogram_multitrial(log(scorelight_var_NI{i}(:, j)), [], [], [], 'Log Principal Position Variance (mm^{2})', ['PM ' num2str(j) ' (Inhibition #' num2str(i) ')']);
            histogram_multitrial(log(speedlight_var_NI{i}(:, j)), [], [], [], 'Log Speed Variance (mm^{2}/s^{2})', ['PM ' num2str(j) ' (Inhibition #' num2str(i) ')']);
            histogram_multitrial(log(accelerationlight_var_NI{i}(:, j)), [], [], [], 'Log Absolute Acceleration Variance (mm^{2}/s^{4})', ['PM ' num2str(j) ' (Inhibition #' num2str(i) ')']);
        end
        
        result.scorelight_mean_NI(i, :) = mean(scorelight_mean_NI{i}, 'omitnan');
        result.speedlight_mean_NI(i, :) = mean(speedlight_mean_NI{i}, 'omitnan');
        result.accelerationlight_mean_NI(i, :) = mean(accelerationlight_mean_NI{i}, 'omitnan');
        
        result.scorelight_var_NI(i, :) = mean(scorelight_var_NI{i}, 'omitnan');
        result.speedlight_var_NI(i, :) = mean(speedlight_var_NI{i}, 'omitnan');
        result.accelerationlight_var_NI(i, :) = mean(accelerationlight_var_NI{i}, 'omitnan');
    end
end

[file, path] = uiputfile('PCAZaxis.mat', 'Save result');
if file ~= 0
    save([path file], 'result');
end