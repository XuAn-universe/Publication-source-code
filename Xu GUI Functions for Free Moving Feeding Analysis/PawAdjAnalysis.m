function PawAdjAnalysis(app, Exp_Path, FrameRate, nmedian)
%*---------------------------------------------------------------------*
% Huang Lab
% Duke University
% Author : Xu An, Jan 2022
% xu.an@duke.edu
% Version: 1.0
%*---------------------------------------------------------------------*
persistent path
curpwd = pwd;
try
    cd(path);
end
path = uigetdir('', 'Pick a folder to save data');
session = split(Exp_Path, '\');
session = session{3};

value = app.TrialsListBox.Value;
value = sort(value);
trials = numel(value);

if app.TrackingDataCheckBox.Value
    label_table = table2array(app.UITable.Data);
else
    helpdlg('''Tracking Data'' needs to be checked');
    return;
end

FrameRate = round(FrameRate);

try
    audiolocation = Exp_Path(1:end-7);
    temp = load([audiolocation '\Detected_Bite_Events.mat']);
    Bite_events = temp.Audio_analysis;
catch
    errordlg('Please detect the bites first.', 'Error');
    return;
end

lcolor = {'r', 'g', 'b', 'c', 'm', 'y', 'k'};
showplot = 1;
for i = 1:trials
    [~, ~, ~, ~, ~, xpawl, ypawl, zpawl, speedpawl, accelerationpawl, ~, ~] = trajectory_postprocessing(10, Exp_Path, value(i), label_table, nmedian, FrameRate);
    [~, ~, ~, ~, ~, xpawr, ypawr, zpawr, speedpawr, accelerationpawr, ~, ~] = trajectory_postprocessing(16, Exp_Path, value(i), label_table, nmedian, FrameRate);
%     [y_bottomPG1, y_bottomPG3, z_bottomPG1, z_bottomPG2, z_bottomPG3, x_bottom, ~, ~, ~, ~, ~, ~] = trajectory_postprocessing(22, Exp_Path, value(i),...
%         label_table, nmedian, FrameRate);
%     [y_centerPG1, y_centerPG3, z_centerPG1, z_centerPG2, z_centerPG3, x_center, ~, ~, ~, ~, ~, ~] = trajectory_postprocessing(33, Exp_Path, value(i),...
%         label_table, nmedian, FrameRate);
    
    frames = size(xpawl, 1);
    t = (1:frames)'/FrameRate;
    
    bite_timestamps = Bite_events(value(i)).time_bites;
    biteIDs = t >= bite_timestamps(1) & t <= bite_timestamps(end);
    
    try
        temp = load([Exp_Path '\LabelledEvents' num2str(value(i)) '.mat']);
        events = temp.LabelledEvents;
        PawLAdjustmentStart = events.PawLAdjustmentStart;
        PawLAdjustmentEnd = events.PawLAdjustmentEnd;
        PawRAdjustmentStart = events.PawRAdjustmentStart;
        PawRAdjustmentEnd = events.PawRAdjustmentEnd;
        ts = min([PawLAdjustmentStart; PawRAdjustmentStart]);
        te = max([PawLAdjustmentEnd; PawRAdjustmentEnd]);
        
        labelspawl = time2category(t, {[PawLAdjustmentStart PawLAdjustmentEnd]}, {'PawLAdjustment'});
        labelspawr = time2category(t, {[PawRAdjustmentStart PawRAdjustmentEnd]}, {'PawRAdjustment'});
    catch
        ts = [];
        te = [];
    end

    pawl = [xpawl(:, 1) ypawl(:, 1) zpawl(:, 1)];
    pawr = [xpawr(:, 1) ypawr(:, 1) zpawr(:, 1)];
    
%     % paw to pasta line distance and paw to pasta bottom distance
%     a =(x_center(:, 1)-x_bottom(:, 1))./(z_centerPG2(:, 1)-z_bottomPG2(:, 1)); % x = a*z+b
%     b = (x_bottom(:, 1).*z_centerPG2(:, 1)-x_center(:, 1).*z_bottomPG2(:, 1))./(z_centerPG2(:, 1)-z_bottomPG2(:, 1));
%     centerPG1 = [a.*z_centerPG1(:, 1)+b y_centerPG1 z_centerPG1];
%     bottomPG1 = [a.*z_bottomPG1(:, 1)+b y_bottomPG1 z_bottomPG1];
%     centerPG3 = [a.*z_centerPG3(:, 1)+b y_centerPG3 z_centerPG3];
%     bottomPG3 = [a.*z_bottomPG3(:, 1)+b y_bottomPG3 z_bottomPG3];
%     pawl2pasta = vecnorm(cross(pawl-centerPG1, pawl-bottomPG1, 2), 2, 2)./vecnorm(centerPG1-bottomPG1, 2, 2);
%     pawl2bottom = dot(pawl-bottomPG1, centerPG1-bottomPG1, 2)./vecnorm(centerPG1-bottomPG1, 2, 2);
%     pawr2pasta = vecnorm(cross(pawr-centerPG3, pawr-bottomPG3, 2), 2, 2)./vecnorm(centerPG3-bottomPG3, 2, 2);
%     pawr2bottom = dot(pawr-bottomPG3, centerPG3-bottomPG3, 2)./vecnorm(centerPG3-bottomPG3, 2, 2);
%     speedpawl2pasta = [nan; diff(pawl2pasta)];
%     speedpawl2bottom = [nan; diff(pawl2bottom)];
%     speedpawr2pasta = [nan; diff(pawr2pasta)];
%     speedpawr2bottom = [nan; diff(pawr2bottom)];
%     pawl2pasta(t <= bite_timestamps(1) | t >= bite_timestamps(end), :) = [];
%     pawl2bottom(t <= bite_timestamps(1) | t >= bite_timestamps(end), :) = [];
%     pawr2pasta(t <= bite_timestamps(1) | t >= bite_timestamps(end), :) = [];
%     pawr2bottom(t <= bite_timestamps(1) | t >= bite_timestamps(end), :) = [];
%     speedpawl2pasta(t <= bite_timestamps(1) | t >= bite_timestamps(end), :) = [];
%     speedpawl2bottom(t <= bite_timestamps(1) | t >= bite_timestamps(end), :) = [];
%     speedpawr2pasta(t <= bite_timestamps(1) | t >= bite_timestamps(end), :) = [];
%     speedpawr2bottom(t <= bite_timestamps(1) | t >= bite_timestamps(end), :) = [];
%     speedpawl2pasta = abs(speedpawl2pasta);
%     speedpawr2pasta = abs(speedpawr2pasta);
    
    pawl2pawr = vecnorm(pawl-pawr, 2, 2);
    speedpawl2pawr = [nan; diff(pawl2pawr)];
    speedpawl = speedpawl(:, 7);
    speedpawr = speedpawr(:, 7);
    accelerationpawl = accelerationpawl(:, 7);
    accelerationpawr = accelerationpawr(:, 7);
    
    % signal normalization
    pawl2pawr = (pawl2pawr-mean(pawl2pawr(biteIDs), 'omitnan'))./std(pawl2pawr(biteIDs), 'omitnan');
    speedpawl2pawr = (speedpawl2pawr-mean(speedpawl2pawr(biteIDs), 'omitnan'))./std(speedpawl2pawr(biteIDs), 'omitnan');
    speedpawl = (speedpawl-mean(speedpawl(biteIDs), 'omitnan'))./std(speedpawl(biteIDs), 'omitnan');
    speedpawr = (speedpawr-mean(speedpawr(biteIDs), 'omitnan'))./std(speedpawr(biteIDs), 'omitnan');
    accelerationpawl = (accelerationpawl-mean(accelerationpawl(biteIDs), 'omitnan'))./std(accelerationpawl(biteIDs), 'omitnan');
    accelerationpawr = (accelerationpawr-mean(accelerationpawr(biteIDs), 'omitnan'))./std(accelerationpawr(biteIDs), 'omitnan');
    
%     [s,w,n] = fsst(fillmissing(pawl2pawr, 'spline'), FrameRate);
    if ~isempty(ts) && ~isempty(te)
        nanID = isnan(speedpawl2pawr); % nan needs to be removed
        
        pawl2pawr(t < ts | t > te | nanID) = [];
        speedpawl2pawr(t < ts | t > te | nanID) = [];
        speedpawl(t < ts | t > te | nanID) = [];
        speedpawr(t < ts | t > te | nanID) = [];
        accelerationpawl(t < ts | t > te | nanID) = [];
        accelerationpawr(t < ts | t > te | nanID) = [];
        
        labelspawl(t < ts | t > te | nanID) = [];
        labelspawr(t < ts | t > te | nanID) = [];
        
        t(t < ts | t > te | nanID) = [];
        
        if ~isempty(labelspawl) && ~isempty(labelspawr)
%             % use this part if there is a problem with class imbalance
%             data2keep = matchlabels(labelspawl, 'PawLAdjustment');
%             data4pawl = {[pawl2pawr(data2keep) speedpawl2pawr(data2keep) speedpawl(data2keep) speedpawr(data2keep)]', labelspawl(data2keep), t(data2keep)'};
%             data2keep = matchlabels(labelspawr, 'PawRAdjustment');
%             data4pawr = {[pawl2pawr(data2keep) speedpawl2pawr(data2keep) speedpawl(data2keep) speedpawr(data2keep)]', labelspawr(data2keep), t(data2keep)'};
            
            data4pawl = {[pawl2pawr speedpawl2pawr speedpawl speedpawr]', labelspawl, t'};
            data4pawr = {[pawl2pawr speedpawl2pawr speedpawl speedpawr]', labelspawr, t'};
            
            save([session ' Trial' num2str(value(i), '%03d') '.mat'], 'data4pawl', 'data4pawr');
        end
    end
    
    if showplot
%         figure;
%         plot(t, pawl2pasta, '-', 'Color', lcolor{1});
%         hold on;
%         plot(t, pawl2bottom, '-', 'Color', lcolor{2});
%         plot(t, pawr2pasta, '-', 'Color', lcolor{3});
%         plot(t, pawr2bottom, '-', 'Color', lcolor{4});
%         legend({'Left paw to pasta', 'Left paw to pasta bottom', 'Right paw to pasta', 'Right paw to pasta bottom'});
%         
%         figure;
%         plot(t, speedpawl2pasta, '-', 'Color', lcolor{1});
%         hold on;
%         plot(t, speedpawl2bottom, '-', 'Color', lcolor{2});
%         plot(t, speedpawr2pasta, '-', 'Color', lcolor{3});
%         plot(t, speedpawr2bottom, '-', 'Color', lcolor{4});
%         legend({'Left paw to pasta speed', 'Left paw to pasta bottom speed', 'Right paw to pasta speed', 'Right paw to pasta bottom speed'});
        
        figure;
        subplot(2, 2, 1);
        plot(t, pawl2pawr, '-', 'Color', lcolor{1});
        set(gca, 'TickDir', 'in', 'FontSize', 12, 'ButtonDownFcn', @extract_figure);
        subplot(2, 2, 2);
        plot(t, speedpawl2pawr, '-', 'Color', lcolor{2});
        set(gca, 'TickDir', 'in', 'FontSize', 12, 'ButtonDownFcn', @extract_figure);
        subplot(2, 2, 3);
        plot(t, speedpawl, '-', 'Color', lcolor{3});
        set(gca, 'TickDir', 'in', 'FontSize', 12, 'ButtonDownFcn', @extract_figure);
        subplot(2, 2, 4);
        plot(t, accelerationpawl, '-', 'Color', lcolor{4});
        set(gca, 'TickDir', 'in', 'FontSize', 12, 'ButtonDownFcn', @extract_figure);
        
        figure;
        plot(t, pawl2pawr, '-', 'Color', lcolor{1});
        hold on;
        plot(t, speedpawl2pawr, '-', 'Color', lcolor{2});
        plot(t, speedpawl, '-', 'Color', lcolor{3});
        plot(t, accelerationpawl, '-', 'Color', lcolor{4});
        legend({'Left paw to right paw', 'Left paw to right paw speed', 'Left paw speed', 'Left paw acceleration'});
        
        figure;
        subplot(2, 2, 1);
        plot(t, pawl2pawr, '-', 'Color', lcolor{1});
        set(gca, 'TickDir', 'in', 'FontSize', 12, 'ButtonDownFcn', @extract_figure);
        subplot(2, 2, 2);
        plot(t, speedpawl2pawr, '-', 'Color', lcolor{2});
        set(gca, 'TickDir', 'in', 'FontSize', 12, 'ButtonDownFcn', @extract_figure);
        subplot(2, 2, 3);
        plot(t, speedpawr, '-', 'Color', lcolor{3});
        set(gca, 'TickDir', 'in', 'FontSize', 12, 'ButtonDownFcn', @extract_figure);
        subplot(2, 2, 4);
        plot(t, accelerationpawr, '-', 'Color', lcolor{4});
        set(gca, 'TickDir', 'in', 'FontSize', 12, 'ButtonDownFcn', @extract_figure);
        
        figure;
        plot(t, pawl2pawr, '-', 'Color', lcolor{1});
        hold on;
        plot(t, speedpawl2pawr, '-', 'Color', lcolor{2});
        plot(t, speedpawr, '-', 'Color', lcolor{3});
        plot(t, accelerationpawr, '-', 'Color', lcolor{4});
        legend({'Left paw to right paw', 'Left paw to right paw speed', 'Right paw speed', 'Right paw acceleration'});
    end
end
cd(curpwd);
msgbox('Done !');