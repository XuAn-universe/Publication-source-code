function [y_PG1, y_PG3, z_PG1, z_PG2, z_PG3, x, y, z, speed, acceleration, laserstart, laserstop] = trajectory_postprocessing(pointID, Exp_Path, trial, label_table, nmedian, FrameRate)
% By Xu An Nov. 2019 @Cold Spring Harbor Laboratory
% xan@cshl.edu
side = 'LR';
if any([1 5 6 7 8 9 10 17] == pointID)
    % left body part
    side = 'L';
end
if any([2 11 12 13 14 15 16 18] == pointID)
    % right body part
    side = 'R';
end
cam1track = find_tracking_filename(Exp_Path, ['PG1camera' num2str(trial)]);
cam2track = find_tracking_filename(Exp_Path, ['PG2camera' num2str(trial)]);
video3location = strsplit(Exp_Path, filesep);
video3location{end-2} = [video3location{end-2} ' Camera3'];
video3location{end-1} = [video3location{end-1} ' PG3'];
video3location = strjoin(video3location, filesep);
cam3track = find_tracking_filename(video3location, ['PG3camera' num2str(trial)]);
cam1track_mm = [cam1track(1:end-4) '.mat'];
cam2track_mm = [cam2track(1:end-4) '.mat'];
cam3track_mm = [cam3track(1:end-4) '.mat'];

load(cam1track_mm, 'num_track');
num_track1_mm = num_track;
load(cam2track_mm, 'num_track');
num_track2_mm = num_track;
load(cam3track_mm, 'num_track');
num_track3_mm = num_track;

y_PG1 = num_track1_mm(:, pointID*3-1);
z_PG1 = num_track1_mm(:, pointID*3);
pcutoff_PG1 = num_track1_mm(:, pointID*3+1);
y_PG1(pcutoff_PG1 <= label_table(pointID, 1)) = NaN;
z_PG1(pcutoff_PG1 <= label_table(pointID, 1)) = NaN;

x_PG2 = num_track2_mm(:, pointID*3-1);
z_PG2 = num_track2_mm(:, pointID*3);
pcutoff_PG2 = num_track2_mm(:, pointID*3+1);
x_PG2(pcutoff_PG2 <= label_table(pointID, 1)) = NaN;
z_PG2(pcutoff_PG2 <= label_table(pointID, 1)) = NaN;

y_PG3 = num_track3_mm(:, pointID*3-1);
z_PG3 = num_track3_mm(:, pointID*3);
pcutoff_PG3 = num_track3_mm(:, pointID*3+1);
y_PG3(pcutoff_PG3 <= label_table(pointID, 1)) = NaN;
z_PG3(pcutoff_PG3 <= label_table(pointID, 1)) = NaN;

% set the location of lost frames to NaN
csvdata1 = load([Exp_Path '\PG1gpio' num2str(trial) '.csv']);
csvdata2 = load([Exp_Path '\PG2gpio' num2str(trial) '.csv']);
csvdata3 = load([video3location '\PG3gpio' num2str(trial) '.csv']);

frameID_PG1 = csvdata1(:, 2)-csvdata1(1, 2)+1;
frameID_PG2 = csvdata2(:, 2)-csvdata2(1, 2)+1;
frameID_PG3 = csvdata3(:, 2)-csvdata3(1, 2)+1;
frames_shared = min([frameID_PG1(end) frameID_PG2(end) frameID_PG3(end)]);

laserstart_PG2 = find(diff([0; csvdata2(:, 1)]) == 1);
laserstart_PG2 = frameID_PG2(laserstart_PG2);
laserstop_PG2 = find(diff([0; csvdata2(:, 1)]) == -1)-1;
laserstop_PG2 = frameID_PG2(laserstop_PG2);
xz_PG2 = NaN(frameID_PG2(end), 2);
xz_PG2(frameID_PG2, 1) = x_PG2;
xz_PG2(frameID_PG2, 2) = z_PG2;
x_PG2 = xz_PG2(1:frames_shared, 1);
z_PG2 = xz_PG2(1:frames_shared, 2);

laserstart_PG1 = find(diff([0; csvdata1(:, 1)]) == 1);
laserstart_PG1 = frameID_PG1(laserstart_PG1);
laserstop_PG1 = find(diff([0; csvdata1(:, 1)]) == -1)-1;
laserstop_PG1 = frameID_PG1(laserstop_PG1);
yz_PG1 = NaN(frameID_PG1(end), 2);
yz_PG1(frameID_PG1, 1) = y_PG1;
yz_PG1(frameID_PG1, 2) = z_PG1;
y_PG1 = yz_PG1(1:frames_shared, 1);
z_PG1 = yz_PG1(1:frames_shared, 2);

laserstart_PG3 = find(diff([0; csvdata3(:, 1)]) == 1);
laserstart_PG3 = frameID_PG3(laserstart_PG3);
laserstop_PG3 = find(diff([0; csvdata3(:, 1)]) == -1)-1;
laserstop_PG3 = frameID_PG3(laserstop_PG3);
yz_PG3 = NaN(frameID_PG3(end), 2);
yz_PG3(frameID_PG3, 1) = y_PG3;
yz_PG3(frameID_PG3, 2) = z_PG3;
y_PG3 = yz_PG3(1:frames_shared, 1);
z_PG3 = yz_PG3(1:frames_shared, 2);

% align the trajectories by linear regression (to correct labeling error)
% z_PG3_as_PG2 = [];
% try
%     mdl = fitlm(z_PG3, z_PG2, 'RobustOpts', 'on');
%     z_PG3_as_PG2 = predict(mdl, z_PG3);
% end
% z_PG1_as_PG2 = [];
% try
%     mdl = fitlm(z_PG1, z_PG2, 'RobustOpts', 'on');
%     z_PG1_as_PG2 = predict(mdl, z_PG1);
% end
% z_PG1_as_PG3 = [];
% try
%     mdl = fitlm(z_PG1, z_PG3, 'RobustOpts', 'on');
%     z_PG1_as_PG3 = predict(mdl, z_PG1);
% end
% z_PG2_as_PG3 = [];
% try
%     mdl = fitlm(z_PG2, z_PG3, 'RobustOpts', 'on');
%     z_PG2_as_PG3 = predict(mdl, z_PG2);
% end
% z_PG2_as_PG1 = [];
% try
%     mdl = fitlm(z_PG2, z_PG1, 'RobustOpts', 'on');
%     z_PG2_as_PG1 = predict(mdl, z_PG2);
% end
% z_PG3_as_PG1 = [];
% try
%     mdl = fitlm(z_PG3, z_PG1, 'RobustOpts', 'on');
%     z_PG3_as_PG1 = predict(mdl, z_PG3);
% end
% y_PG1_as_PG3 = [];
% try
%     mdl = fitlm(y_PG1, y_PG3, 'RobustOpts', 'on');
%     y_PG1_as_PG3 = predict(mdl, y_PG1);
% end
% y_PG3_as_PG1 = [];
% try
%     mdl = fitlm(y_PG3, y_PG1, 'RobustOpts', 'on');
%     y_PG3_as_PG1 = predict(mdl, y_PG3);
% end

xavg = x_PG2;
switch side
    case 'L'
        yavg = y_PG1;
        yavg(isnan(yavg)) = y_PG3(isnan(yavg));
        zavg = z_PG1;
        z_PG2(isnan(z_PG2)) = z_PG3(isnan(z_PG2));
        zavg(isnan(zavg)) = z_PG2(isnan(zavg));
    case 'R'
        yavg = y_PG3;
        yavg(isnan(yavg)) = y_PG1(isnan(yavg));
        zavg = z_PG3;
        z_PG2(isnan(z_PG2)) = z_PG1(isnan(z_PG2));
        zavg(isnan(zavg)) = z_PG2(isnan(zavg));
    case 'LR'
        yavg = mean([y_PG1 y_PG3], 2, 'omitnan');
        zavg = z_PG2;
        z_PG1PG3 = mean([z_PG1 z_PG3], 2, 'omitnan');
        zavg(isnan(zavg)) = z_PG1PG3(isnan(zavg));
end

% filter the trajectories with median filter
xavg_filtered = medfilt1(xavg, nmedian, 'omitnan', 'truncate');
yavg_filtered = medfilt1(yavg, nmedian, 'omitnan', 'truncate');
zavg_filtered = medfilt1(zavg, nmedian, 'omitnan', 'truncate');

% speed & acceleration
speedx = diff(xavg)*FrameRate;
speedx = [NaN; speedx];
accelerationx = diff(abs(speedx))*FrameRate;
accelerationx = [NaN; accelerationx];
% speedx = medfilt1(speedx, nmedian, 'omitnan', 'truncate');
% accelerationx = medfilt1(accelerationx, nmedian, 'omitnan', 'truncate');

speedy = diff(yavg)*FrameRate;
speedy = [NaN; speedy];
accelerationy = diff(abs(speedy))*FrameRate;
accelerationy = [NaN; accelerationy];
% speedy = medfilt1(speedy, nmedian, 'omitnan', 'truncate');
% accelerationy = medfilt1(accelerationy, nmedian, 'omitnan', 'truncate');

speedz = diff(zavg)*FrameRate;
speedz = [NaN; speedz];
accelerationz = diff(abs(speedz))*FrameRate;
accelerationz = [NaN; accelerationz];
% speedz = medfilt1(speedz, nmedian, 'omitnan', 'truncate');
% accelerationz = medfilt1(accelerationz, nmedian, 'omitnan', 'truncate');

speedxy = sqrt(diff(xavg).^2+diff(yavg).^2)*FrameRate;
speedxy = [NaN; speedxy];
accelerationxy = diff(speedxy)*FrameRate;
accelerationxy = [NaN; accelerationxy];
% speedxy = medfilt1(speedxy, nmedian, 'omitnan', 'truncate');
% accelerationxy = medfilt1(accelerationxy, nmedian, 'omitnan', 'truncate');

speedxz = sqrt(diff(xavg).^2+diff(zavg).^2)*FrameRate;
speedxz = [NaN; speedxz];
accelerationxz = diff(speedxz)*FrameRate;
accelerationxz = [NaN; accelerationxz];
% speedxz = medfilt1(speedxz, nmedian, 'omitnan', 'truncate');
% accelerationxz = medfilt1(accelerationxz, nmedian, 'omitnan', 'truncate');

speedyz = sqrt(diff(yavg).^2+diff(zavg).^2)*FrameRate;
speedyz = [NaN; speedyz];
accelerationyz = diff(speedyz)*FrameRate;
accelerationyz = [NaN; accelerationyz];
% speedyz = medfilt1(speedyz, nmedian, 'omitnan', 'truncate');
% accelerationyz = medfilt1(accelerationyz, nmedian, 'omitnan', 'truncate');

speedxyz = sqrt(diff(xavg).^2+diff(yavg).^2+diff(zavg).^2)*FrameRate;
speedxyz = [NaN; speedxyz];
accelerationxyz = diff(speedxyz)*FrameRate;
accelerationxyz = [NaN; accelerationxyz];
% speedxyz = medfilt1(speedxyz, nmedian, 'omitnan', 'truncate');
% accelerationxyz = medfilt1(accelerationxyz, nmedian, 'omitnan', 'truncate');

laserstart = min([laserstart_PG1 laserstart_PG2 laserstart_PG3], [], 2);
laserstop = max([laserstop_PG1 laserstop_PG2 laserstop_PG3], [], 2);
if numel(laserstart)-numel(laserstop) == 1
    laserstart(end) = [];
end

x = [xavg xavg_filtered];
y = [yavg yavg_filtered];
z = [zavg zavg_filtered];
speed = [speedx speedy speedz speedxy speedxz speedyz speedxyz];
acceleration = [accelerationx accelerationy accelerationz accelerationxy accelerationxz accelerationyz accelerationxyz];

    function filename = find_tracking_filename(Exp_Path, camID)
        filename = [];
        length_camID = length(camID);
        curpwd = pwd;
        cd(Exp_Path);
        DirList = dir;
        for i = 3:size(DirList)
            if numel(DirList(i).name) > length_camID
                if strcmp(DirList(i).name(1:length_camID), camID) && strcmp(DirList(i).name(end-3:end), '.csv') && ~isempty(strfind(DirList(i).name, 'DeepCut'))
                    filename = [Exp_Path '\' DirList(i).name];
                end
            end
        end
        cd(curpwd);
    end
end