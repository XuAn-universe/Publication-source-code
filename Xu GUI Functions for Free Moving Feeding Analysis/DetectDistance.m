function DetectDistance(app, Exp_Path, FrameRate, nmedian)
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

detectdist = [];
detectdist_all = cell(1, numel(Exp_Paths));

FrameRate = round(FrameRate);

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
    for j = 1:numel(Video_annotation)
        detectdist_all{i}(j).data = [];
        if ~Video_annotation(j).Disgard
            [y_topPG1, y_topPG3, z_topPG1, z_topPG2, z_topPG3, x_top, ~, ~, ~, ~, laserstart, laserstop] = trajectory_postprocessing(21, Exp_Paths{i}, j,...
                label_table, nmedian, FrameRate);
            [y_bottomPG1, y_bottomPG3, z_bottomPG1, z_bottomPG2, z_bottomPG3, x_bottom, ~, ~, ~, ~, ~, ~] = trajectory_postprocessing(22, Exp_Paths{i}, j,...
                label_table, nmedian, FrameRate);
            [y_centerPG1, y_centerPG3, z_centerPG1, z_centerPG2, z_centerPG3, x_center, ~, ~, ~, ~, ~, ~] = trajectory_postprocessing(33, Exp_Paths{i}, j,...
                label_table, nmedian, FrameRate);
            [~, ~, ~, ~, ~, x_nose, y_nose, z_nose, ~, ~, ~, ~] = trajectory_postprocessing(3, Exp_Paths{i}, j,...
                label_table, nmedian, FrameRate);
            
            timestamps = [];
            try
                temp = load([Exp_Paths{i} '\LabelledEvents' num2str(j) '.mat']);
                LabelledEvents = temp.LabelledEvents;
                timestamps = union(LabelledEvents.RetrievalStart, LabelledEvents.MouthRetrievalStart);
            end
            if isempty(timestamps)
                continue;
            else
                eventframe = round(timestamps*FrameRate);
            end
            
            xyz_topPG1 = [x_top(eventframe, 1) y_topPG1(eventframe, 1) mean([z_topPG1(eventframe, 1) z_topPG2(eventframe, 1)], 2)];
            xyz_bottomPG1 = [x_bottom(eventframe, 1) y_bottomPG1(eventframe, 1) mean([z_bottomPG1(eventframe, 1) z_bottomPG2(eventframe, 1)], 2)];
            xyz_centerPG1 = [x_center(eventframe, 1) y_centerPG1(eventframe, 1) mean([z_centerPG1(eventframe, 1) z_centerPG2(eventframe, 1)], 2)];
            xyz_topPG3 = [x_top(eventframe, 1) y_topPG3(eventframe, 1) mean([z_topPG3(eventframe, 1) z_topPG2(eventframe, 1)], 2)];
            xyz_bottomPG3 = [x_bottom(eventframe, 1) y_bottomPG3(eventframe, 1) mean([z_bottomPG3(eventframe, 1) z_bottomPG2(eventframe, 1)], 2)];
            xyz_centerPG3 = [x_center(eventframe, 1) y_centerPG3(eventframe, 1) mean([z_centerPG3(eventframe, 1) z_centerPG2(eventframe, 1)], 2)];
            xyz_nose = [x_nose(eventframe, 1) y_nose(eventframe, 1) z_nose(eventframe, 1)];
            distance = zeros(size(xyz_nose, 1), 6);
            distance(:, 1) = sqrt(sum((xyz_nose-xyz_topPG1).^2, 2));
            distance(:, 2) = sqrt(sum((xyz_nose-xyz_bottomPG1).^2, 2));
            distance(:, 3) = sqrt(sum((xyz_nose-xyz_centerPG1).^2, 2));
            distance(:, 4) = sqrt(sum((xyz_nose-xyz_topPG3).^2, 2));
            distance(:, 5) = sqrt(sum((xyz_nose-xyz_bottomPG3).^2, 2));
            distance(:, 6) = sqrt(sum((xyz_nose-xyz_centerPG3).^2, 2));
            for k = 1:size(distance, 1)
                if any(~isnan(distance(k, :)))
                    [minValue, minID] = min(distance(k, :));
                    if minValue < 10
                        detectdist = [detectdist minValue];
                        detectdist_all{i}(j).data = [detectdist_all{i}(j).data minValue];
                        im = retrieveframe(app, Exp_Paths{i}, j, timestamps(k));
                        switch minID
                            case 1
                                im(:, :, :, 2) = insertText(im(:, :, :, 2), [0 0], ['Nose to Top PG1 = ' num2str(minValue, '%.1f')], 'TextColor', 'white', 'BoxColor', 'black', 'FontSize', 30);
                            case 2
                                im(:, :, :, 2) = insertText(im(:, :, :, 2), [0 0], ['Nose to Bottom PG1 = ' num2str(minValue, '%.1f')], 'TextColor', 'white', 'BoxColor', 'black', 'FontSize', 30);
                            case 3
                                im(:, :, :, 2) = insertText(im(:, :, :, 2), [0 0], ['Nose to Center PG1 = ' num2str(minValue, '%.1f')], 'TextColor', 'white', 'BoxColor', 'black', 'FontSize', 30);
                            case 4
                                im(:, :, :, 2) = insertText(im(:, :, :, 2), [0 0], ['Nose to Top PG3 = ' num2str(minValue, '%.1f')], 'TextColor', 'white', 'BoxColor', 'black', 'FontSize', 30);
                            case 5
                                im(:, :, :, 2) = insertText(im(:, :, :, 2), [0 0], ['Nose to Bottom PG3 = ' num2str(minValue, '%.1f')], 'TextColor', 'white', 'BoxColor', 'black', 'FontSize', 30);
                            case 6
                                im(:, :, :, 2) = insertText(im(:, :, :, 2), [0 0], ['Nose to Center PG3 = ' num2str(minValue, '%.1f')], 'TextColor', 'white', 'BoxColor', 'black', 'FontSize', 30);
                        end
                        
                        if ~writetiff
                            imwrite([im(:, :, :, 1) im(:, :, :, 2) im(:, :, :, 3)], 'ValidationVideo.tiff', 'tiff', 'Compression', 'none');
                            writetiff = true;
                        else
                            imwrite([im(:, :, :, 1) im(:, :, :, 2) im(:, :, :, 3)], 'ValidationVideo.tiff', 'tiff', 'Compression', 'none', 'WriteMode', 'append');
                        end
                    end
                end
            end
        end
        workbar(j/numel(Video_annotation), [num2str(j) '/' num2str(numel(Video_annotation))], 'Progress'); 
    end
end

assignin('base', 'detectdist_all', detectdist_all);
assignin('base', 'detectdist', detectdist);
save('detectdistance.mat', 'detectdist');
cd(curpwd);
msgbox('Done !');
end