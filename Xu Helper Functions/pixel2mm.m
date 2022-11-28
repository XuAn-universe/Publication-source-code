function pixel2mm(Exp_Path, CP_Path)
% eg. pixel2mm('F:\Free-moving Feeding Data\072119_AngelHair_Fezf2CreER_GtACR1-8_Inhibition\exp001', 'F:\Camera Parameters');
% eg. pixel2mm('F:\Free-moving Feeding Data\072119_AngelHair_Fezf2CreER_GtACR1-8_Inhibition\exp001', 'F:\Camera Parameters @duke');
% By Xu An Nov. 2019 @Cold Spring Harbor Laboratory
% xan@cshl.edu
tic;
inch2mmcst = 25.4;
load([CP_Path '\ref_point.mat'], 'ref_point');
ref_points = ref_point;
load([CP_Path '\cameraParams_PG1.mat'], 'cameraParams');
CPs{1} = cameraParams;
load([CP_Path '\cameraParams_PG2.mat'], 'cameraParams');
CPs{2} = cameraParams;
load([CP_Path '\cameraParams_PG3.mat'], 'cameraParams');
CPs{3} = cameraParams;

curpwd = pwd;
for i = 1:3
    if i ~= 3
        CSV_Path = Exp_Path;
    else
        CSV_Path = strsplit(Exp_Path, filesep);
        CSV_Path{end-2} = [CSV_Path{end-2} ' Camera3'];
        CSV_Path{end-1} = [CSV_Path{end-1} ' PG3'];
        CSV_Path = strjoin(CSV_Path, filesep);
    end
    cd(CSV_Path);
    DirList = dir;
    n = 0;
    filename = cell(0);
    camID = ['PG' num2str(i) 'camera'];
    for j = 3:size(DirList)
        if numel(DirList(j).name) > 20
            if ~isempty(strfind(DirList(j).name, camID)) && strcmp(DirList(j).name(end-3:end), '.csv') && ~isempty(strfind(DirList(j).name, 'DeepCut_resnet50'))
                n = n+1;
                filename{n} = [CSV_Path '\' DirList(j).name];
            end
        end
    end
    
    CP = CPs{i};
    ref = ref_points(i, :);
    ref = undistortPoints(ref, CP);
    ref = pointsToWorld(CP, CP.RotationMatrices(:, :, end), CP.TranslationVectors(end, :), ref);
    for j = 1:n
        [num_track, ~, ~] = xlsread(filename{j});
        points = (size(num_track, 2)-1)/3;
        temp = cell(1, points);
        for k = 1:points
            temp{k} = num_track(:, k*3-1:k*3);
        end
        parfor k = 1:points
            frames = size(temp{k}, 1);
            nn = 0;
            while frames > 0
                if frames >= 50
                    temp{k}(1+nn*50:(nn+1)*50, :) = undistortPoints(temp{k}(1+nn*50:(nn+1)*50, :), CP);
                else
                    temp{k}(1+nn*50:end, :) = undistortPoints(temp{k}(1+nn*50:end, :), CP);
                end
                frames = frames-50;
                nn = nn+1;
            end
            temp{k} = pointsToWorld(CP, CP.RotationMatrices(:, :, end), CP.TranslationVectors(end, :), temp{k});
            temp{k}(:, 1) = temp{k}(:, 1)-ref(1);
            temp{k}(:, 2) = temp{k}(:, 2)-ref(2);
            if i ~= 3
                temp{k}(:, 2) = -temp{k}(:, 2);
            end
        end
        for k = 1:points
            num_track(:, k*3-1:k*3) = temp{k}*inch2mmcst;
        end
        save([filename{j}(1:end-4) '.mat'], 'num_track');
        disp(filename{j}(1:end-4));
    end
end
cd(curpwd);
msgbox('Done !');
toc