function timedt = FoodDetection(app, Exp_Path, FrameRate, nmedian)
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

timedt_NI = [];
timedt_I = [];
% timefb_NI = [];
% timefb_I = [];
timedt_all = cell(1, numel(Exp_Paths));
timedt = [];

FrameRate = round(FrameRate);
fooddetect_threshold = app.thrEditField.Value;

curpwd = pwd;
SaveDir = uigetdir(Exp_Path, 'Pick a folder to save images for validation');
if ~SaveDir
    return;
else
    cd(SaveDir);
end
writetiff = false;

for i = 1:numel(Exp_Paths)
    load([Exp_Paths{i} '\Analysis_Session.mat'], 'Video_annotation');
    workbar(0, 'Computing Ongoing...', 'Progress'); 
%     try
%         audiolocation = Exp_Paths{i}(1:end-7);
%         temp = load([audiolocation '\Detected_Bite_Events.mat']);
%         Bite_events = temp.Audio_analysis;
%     catch
%         errordlg('Please detect the bites first.', 'Error');
%         return;
%     end
    for j = 1:numel(Video_annotation)
        timedt_all{i}(j) = NaN;
        if ~Video_annotation(j).Disgard
            [y_topPG1, y_topPG3, z_topPG1, z_topPG2, z_topPG3, x_top, ~, ~, ~, ~, laserstart, laserstop] = trajectory_postprocessing(21, Exp_Paths{i}, j,...
                label_table, nmedian, FrameRate);
            [y_bottomPG1, y_bottomPG3, z_bottomPG1, z_bottomPG2, z_bottomPG3, x_bottom, ~, ~, ~, ~, ~, ~] = trajectory_postprocessing(22, Exp_Paths{i}, j,...
                label_table, nmedian, FrameRate);
            [y_centerPG1, y_centerPG3, z_centerPG1, z_centerPG2, z_centerPG3, x_center, ~, ~, ~, ~, ~, ~] = trajectory_postprocessing(33, Exp_Paths{i}, j,...
                label_table, nmedian, FrameRate);
            [~, ~, ~, ~, ~, x_nose, y_nose, z_nose, ~, ~, ~, ~] = trajectory_postprocessing(3, Exp_Paths{i}, j,...
                label_table, nmedian, FrameRate);
            [~, ~, ~, ~, ~, x_leftpaw, y_leftpaw, z_leftpaw, ~, ~, ~, ~] = trajectory_postprocessing(10, Exp_Paths{i}, j,...
                label_table, nmedian, FrameRate);
            
            xyz_topPG1 = [x_top(:, 1) y_topPG1(:, 1) mean([z_topPG1(:, 1) z_topPG2(:, 1)], 2)];
            xyz_bottomPG1 = [x_bottom(:, 1) y_bottomPG1(:, 1) mean([z_bottomPG1(:, 1) z_bottomPG2(:, 1)], 2)];
            xyz_centerPG1 = [x_center(:, 1) y_centerPG1(:, 1) mean([z_centerPG1(:, 1) z_centerPG2(:, 1)], 2)];
            xyz_topPG3 = [x_top(:, 1) y_topPG3(:, 1) mean([z_topPG3(:, 1) z_topPG2(:, 1)], 2)];
            xyz_bottomPG3 = [x_bottom(:, 1) y_bottomPG3(:, 1) mean([z_bottomPG3(:, 1) z_bottomPG2(:, 1)], 2)];
            xyz_centerPG3 = [x_center(:, 1) y_centerPG3(:, 1) mean([z_centerPG3(:, 1) z_centerPG2(:, 1)], 2)];
            
            frame_leftpaw = find(~isnan(x_leftpaw(:, 1)) & ~isnan(y_leftpaw(:, 1)) & ~isnan(z_leftpaw(:, 1)), 1, 'first');
            visible = 0;
            for k = 1:frame_leftpaw
                if all(~isnan(xyz_topPG1(k, :))) || all(~isnan(xyz_bottomPG1(k, :))) || all(~isnan(xyz_centerPG1(k, :))) || all(~isnan(xyz_topPG3(k, :))) || all(~isnan(xyz_bottomPG3(k, :))) || all(~isnan(xyz_centerPG3(k, :)))
                    visible = 1;
                    break;
                end
            end
            if visible
                xyz_nose = [x_nose(:, 1) y_nose(:, 1) z_nose(:, 1)];
                distance = zeros(size(xyz_nose, 1), 6);
                distance(:, 1) = sqrt(sum((xyz_nose-xyz_topPG1).^2, 2));
                distance(:, 2) = sqrt(sum((xyz_nose-xyz_bottomPG1).^2, 2));
                distance(:, 3) = sqrt(sum((xyz_nose-xyz_centerPG1).^2, 2));
                distance(:, 4) = sqrt(sum((xyz_nose-xyz_topPG3).^2, 2));
                distance(:, 5) = sqrt(sum((xyz_nose-xyz_bottomPG3).^2, 2));
                distance(:, 6) = sqrt(sum((xyz_nose-xyz_centerPG3).^2, 2));
                frame_detect = find(sum(distance <= fooddetect_threshold, 2), 1, 'first');
                if isempty(frame_detect)
                    continue;
                end
                disptext = [];
                for k = find(distance(frame_detect, :) <= fooddetect_threshold)
                    switch k
                        case 1
                            disptext = [disptext 'PG1 top; '];
                        case 2
                            disptext = [disptext 'PG1 bottom; '];
                        case 3
                            disptext = [disptext 'PG1 center; '];
                        case 4
                            disptext = [disptext 'PG3 top; '];
                        case 5
                            disptext = [disptext 'PG3 bottom; '];
                        case 6
                            disptext = [disptext 'PG3 center; '];
                    end
                end
                time_detect = frame_detect/FrameRate;
                disptext = [disptext num2str(time_detect, '%.2f') ' s'];
                im = retrieveframe(app, Exp_Paths{i}, j, time_detect);
                im(:, :, :, 2) = insertText(im(:, :, :, 2), [0 0], disptext, 'TextColor', 'white', 'BoxColor', 'black', 'FontSize', 30);
                if ~writetiff
                    imwrite([im(:, :, :, 1) im(:, :, :, 2) im(:, :, :, 3)], 'FoodDetection.tiff', 'tiff', 'Compression', 'none');
                    writetiff = true;
                else
                    imwrite([im(:, :, :, 1) im(:, :, :, 2) im(:, :, :, 3)], 'FoodDetection.tiff', 'tiff', 'Compression', 'none', 'WriteMode', 'append');
                end
                
                timedt_all{i}(j) = time_detect;
                
%                 bite_timestamps = Bite_events(j).time_bites;
%                 time_firstbite = bite_timestamps(1);
                if app.NoLightButton.Value
                    laserstart = [];
                    laserstop = [];
                else
%                     laser_timestamps = Bite_events(j).laser_timestamps;
                end
                if ~isempty(laserstart)
                    timedt_I(end+1) = time_detect;
                    timedt = [timedt; time_detect 1];
%                     timefb_I(end+1) = time_firstbite;
                else
                    timedt_NI(end+1) = time_detect;
                    timedt = [timedt; time_detect 0];
%                     timefb_NI(end+1) = time_firstbite;
                end
            end
        end
        workbar(j/numel(Video_annotation), [num2str(j) '/' num2str(numel(Video_annotation))], 'Progress'); 
    end
end
assignin('base', 'timedt_all', timedt_all);
assignin('base', 'timedt', timedt);
save('timedt.mat', 'timedt');
cd(curpwd);
msgbox('Done !');

if ~app.NoLightButton.Value
    figure;
    hold on;
    bar(1, mean(timedt_NI), 0.8, 'FaceColor', [0.5 0.5 0.5]);
    plot(1+rand(1, numel(timedt_NI))*0.8-0.4, timedt_NI, 'o',...
        'MarkerFaceColor', [0 0 0], 'MarkerEdgeColor', 'none');
    bar(2, mean(timedt_I), 0.8, 'FaceColor', [0 1 0]);
    plot(2+rand(1, numel(timedt_I))*0.8-0.4, timedt_I, 'o',...
        'MarkerFaceColor', [0 0.75 0], 'MarkerEdgeColor', 'none');
    errorbar([1 2], [mean(timedt_NI) mean(timedt_I)],...
        [std(timedt_NI)/sqrt(numel(timedt_NI)) std(timedt_I)/sqrt(numel(timedt_I))],...
        'k', 'LineStyle', 'none');
    set(gca, 'XTick', [1 2], 'XTickLabel', {'Control', 'Inhibition'}, 'XLim', [0.4 2.6], 'TickLength', [0 0], 'FontSize', 12);
    ylabel('Detection Time (s)');
    fprintf('\n');
    p = ranksum(timedt_NI, timedt_I, 'tail', 'both');
    disp(['detection time (both): p = ' num2str(p)]);
    p = ranksum(timedt_NI, timedt_I, 'tail', 'left');
    disp(['detection time (left): p = ' num2str(p)]);
end

% timepk_NI = timefb_NI-timedt_NI;
% timepk_I = timefb_I-timedt_I;
% subplot(1, 2, 2);
% hold on;
% bar(1, mean(timepk_NI), 0.8, 'FaceColor', [0.5 0.5 0.5]);
% plot(1+rand(1, numel(timepk_NI))*0.8-0.4, timepk_NI, 'o',...
%     'MarkerFaceColor', [0 0 0], 'MarkerEdgeColor', 'none');
% bar(2, mean(timepk_I), 0.8, 'FaceColor', [0 1 0]);
% plot(2+rand(1, numel(timepk_I))*0.8-0.4, timepk_I, 'o',...
%     'MarkerFaceColor', [0 0.75 0], 'MarkerEdgeColor', 'none');
% errorbar([1 2], [mean(timepk_NI) mean(timepk_I)],...
%     [std(timepk_NI)/sqrt(numel(timepk_NI)) std(timepk_I)/sqrt(numel(timepk_I))],...
%     'k', 'LineStyle', 'none');
% set(gca, 'XTick', [1 2], 'XTickLabel', {'Control', 'Inhibition'}, 'XLim', [0.4 2.6], 'TickLength', [0 0], 'FontSize', 12);
% ylabel('Picking Up Time (s)');
% fprintf('\n');
% p = ranksum(timepk_NI, timepk_I, 'tail', 'both');
% disp(['picking up time (both): p = ' num2str(p)]);
% p = ranksum(timepk_NI, timepk_I, 'tail', 'left');
% disp(['picking up time (left): p = ' num2str(p)]);