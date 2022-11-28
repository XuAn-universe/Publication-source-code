%*---------------------------------------------------------------------*
% Huang Lab
% Duke University
% Author : Xu An, May 2022
% xu.an@duke.edu
% Version: 1.0
%*---------------------------------------------------------------------*
%% hand adjustment posture analysis across mice
curpwd = pwd;
try
   cd(pathname); 
end
[filename, pathname] = uigetfile('*.mat', 'Pick all of the data set', 'MultiSelect', 'on');
if isequal(filename, 0)
    cd(curpwd);
    return;
end

N = numel(filename);

clc;
FrameRate = 120;
trange = 0.5;
trange = trange*FrameRate;

thr = 15;
pcheck = 1;
stp = 2;
edges = 0:stp:180;

% cnose_NI = zeros(3, N, 3);
cnose_NI = nan(3, N, 3);
ceyel_NI = nan(3, N, 3);
% ceyer_NI = nan(3, N, 3);
ceyer_NI = zeros(3, N, 3);
cpawl_NI = nan(3, N, 3);
cpawr_NI = nan(3, N, 3);
ctop_NI = nan(3, N, 3);
cbottom_NI = nan(3, N, 3);
side_NI = nan(N, 2);
orientation_all_NI.xy = [];
orientation_all_NI.xz = [];
orientation_all_NI.yz = [];
orientation_all_NI.xy_avg = nan(1, N);
orientation_all_NI.xz_avg = nan(1, N);
orientation_all_NI.yz_avg = nan(1, N);
pawgrasp_L_NI = nan(3, N);
pawguide_L_NI = nan(3, N);
pawgrasp_R_NI = nan(3, N);
pawguide_R_NI = nan(3, N);

% cnose_I = zeros(3, N, 3);
cnose_I = nan(3, N, 3);
ceyel_I = nan(3, N, 3);
% ceyer_I = nan(3, N, 3);
ceyer_I = zeros(3, N, 3);
cpawl_I = nan(3, N, 3);
cpawr_I = nan(3, N, 3);
ctop_I = nan(3, N, 3);
cbottom_I = nan(3, N, 3);
side_I = nan(N, 2);
orientation_all_I.xy = [];
orientation_all_I.xz = [];
orientation_all_I.yz = [];
orientation_all_I.xy_avg = nan(1, N);
orientation_all_I.xz_avg = nan(1, N);
orientation_all_I.yz_avg = nan(1, N);
pawgrasp_L_I = nan(3, N);
pawguide_L_I = nan(3, N);
pawgrasp_R_I = nan(3, N);
pawguide_R_I = nan(3, N);

xmidline = nan(1, N);
corrX_pawr.rsquare = nan(N, 4);
corrX_pawr.slope = nan(N, 4);
corrX_pawr.rho = nan(N, 4);
corrX_pawr.p = nan(N, 4);
corrX_pawl.rsquare = nan(N, 4);
corrX_pawl.slope = nan(N, 4);
corrX_pawl.rho = nan(N, 4);
corrX_pawl.p = nan(N, 4);
corrY_pawr.rsquare = nan(N, 4);
corrY_pawr.slope = nan(N, 4);
corrY_pawr.rho = nan(N, 4);
corrY_pawr.p = nan(N, 4);
corrY_pawl.rsquare = nan(N, 4);
corrY_pawl.slope = nan(N, 4);
corrY_pawl.rho = nan(N, 4);
corrY_pawl.p = nan(N, 4);
corrZ_pawr.rsquare = nan(N, 4);
corrZ_pawr.slope = nan(N, 4);
corrZ_pawr.rho = nan(N, 4);
corrZ_pawr.p = nan(N, 4);
corrZ_pawl.rsquare = nan(N, 4);
corrZ_pawl.slope = nan(N, 4);
corrZ_pawl.rho = nan(N, 4);
corrZ_pawl.p = nan(N, 4);

for i = 1:N
    temp = load([pathname filename{i}]);
    result = temp.result;
    if ~(sum(result.IDleft_front_NI-result.IDleft_top_NI) == 0 && sum(result.IDright_front_NI-result.IDright_top_NI) == 0)
        errordlg(['Check ' filename{i} ' for the left and right IDs']);
        return;
    end
    IDleft = result.IDleft_front_NI;
    IDright = result.IDright_front_NI;
    IDall = [IDleft IDright];
    side_NI(i, 1) = numel(IDleft);
    side_NI(i, 2) = numel(IDright);
    
    ctopNI = squeeze(result.ctop_all_NI(trange+1, :, IDall));
    cbottomNI = squeeze(result.cbottom_all_NI(trange+1, :, IDall));
    cpawlNI = squeeze(result.cpawl_all_NI(trange+1, :, IDall));
    cpawrNI = squeeze(result.cpawr_all_NI(trange+1, :, IDall));
    ceyelNI = squeeze(result.ceyel_all_NI(trange+1, 1, IDall));
    
    if numel(IDall) > thr
        ceyel_NI(:, i, 1) = mean(squeeze(result.ceyel_all_NI(trange+1, :, IDall)), 2, 'omitnan');
%         ceyer_NI(:, i, 1) = mean(squeeze(result.ceyer_all_NI(trange+1, :, IDall)), 2, 'omitnan');
        cnose_NI(:, i, 1) = mean(squeeze(result.cnose_all_NI(trange+1, :, IDall)), 2, 'omitnan');
        cpawl_NI(:, i, 1) = mean(squeeze(result.cpawl_all_NI(trange+1, :, IDall)), 2, 'omitnan');
        cpawr_NI(:, i, 1) = mean(squeeze(result.cpawr_all_NI(trange+1, :, IDall)), 2, 'omitnan');
        ctop_NI(:, i, 1) = mean(squeeze(result.ctop_all_NI(trange+1, :, IDall)), 2, 'omitnan');
        cbottom_NI(:, i, 1) = mean(squeeze(result.cbottom_all_NI(trange+1, :, IDall)), 2, 'omitnan');
    end
    if numel(IDleft) > thr
        ceyel_NI(:, i, 2) = mean(squeeze(result.ceyel_all_NI(trange+1, :, IDleft)), 2, 'omitnan');
%         ceyer_NI(:, i, 2) = mean(squeeze(result.ceyer_all_NI(trange+1, :, IDleft)), 2, 'omitnan');
        cnose_NI(:, i, 2) = mean(squeeze(result.cnose_all_NI(trange+1, :, IDleft)), 2, 'omitnan');
        cpawl_NI(:, i, 2) = mean(squeeze(result.cpawl_all_NI(trange+1, :, IDleft)), 2, 'omitnan');
        cpawr_NI(:, i, 2) = mean(squeeze(result.cpawr_all_NI(trange+1, :, IDleft)), 2, 'omitnan');
        ctop_NI(:, i, 2) = mean(squeeze(result.ctop_all_NI(trange+1, :, IDleft)), 2, 'omitnan');
        cbottom_NI(:, i, 2) = mean(squeeze(result.cbottom_all_NI(trange+1, :, IDleft)), 2, 'omitnan');
        
%         pawgrasp_L_NI(1, i) = mean(abs(result.cpawl_all_NI(trange+1, 1, IDleft)-result.ceyel_all_NI(trange+1, 1, IDleft)/2), 'omitnan');
        pawgrasp_L_NI(1, i) = mean(result.cpawl_all_NI(trange+1, 1, IDleft), 'omitnan');
        pawgrasp_L_NI(2, i) = mean(result.cpawl_all_NI(trange+1, 2, IDleft), 'omitnan');
        pawgrasp_L_NI(3, i) = mean(result.cpawl_all_NI(trange+1, 3, IDleft), 'omitnan');
%         pawguide_L_NI(1, i) = mean(abs(result.cpawr_all_NI(trange+1, 1, IDleft)-result.ceyel_all_NI(trange+1, 1, IDleft)/2), 'omitnan');
        pawguide_L_NI(1, i) = mean(result.cpawr_all_NI(trange+1, 1, IDleft), 'omitnan');
        pawguide_L_NI(2, i) = mean(result.cpawr_all_NI(trange+1, 2, IDleft), 'omitnan');
        pawguide_L_NI(3, i) = mean(result.cpawr_all_NI(trange+1, 3, IDleft), 'omitnan');
    end
    if numel(IDright) > thr
        ceyel_NI(:, i, 3) = mean(squeeze(result.ceyel_all_NI(trange+1, :, IDright)), 2, 'omitnan');
%         ceyer_NI(:, i, 3) = mean(squeeze(result.ceyer_all_NI(trange+1, :, IDright)), 2, 'omitnan');
        cnose_NI(:, i, 3) = mean(squeeze(result.cnose_all_NI(trange+1, :, IDright)), 2, 'omitnan');
        cpawl_NI(:, i, 3) = mean(squeeze(result.cpawl_all_NI(trange+1, :, IDright)), 2, 'omitnan');
        cpawr_NI(:, i, 3) = mean(squeeze(result.cpawr_all_NI(trange+1, :, IDright)), 2, 'omitnan');
        ctop_NI(:, i, 3) = mean(squeeze(result.ctop_all_NI(trange+1, :, IDright)), 2, 'omitnan');
        cbottom_NI(:, i, 3) = mean(squeeze(result.cbottom_all_NI(trange+1, :, IDright)), 2, 'omitnan');
        
%         pawgrasp_R_NI(1, i) = mean(abs(result.cpawr_all_NI(trange+1, 1, IDright)-result.ceyel_all_NI(trange+1, 1, IDright)/2), 'omitnan');
        pawgrasp_R_NI(1, i) = mean(result.cpawr_all_NI(trange+1, 1, IDright), 'omitnan');
        pawgrasp_R_NI(2, i) = mean(result.cpawr_all_NI(trange+1, 2, IDright), 'omitnan');
        pawgrasp_R_NI(3, i) = mean(result.cpawr_all_NI(trange+1, 3, IDright), 'omitnan');
%         pawguide_R_NI(1, i) = mean(abs(result.cpawl_all_NI(trange+1, 1, IDright)-result.ceyel_all_NI(trange+1, 1, IDright)/2), 'omitnan');
        pawguide_R_NI(1, i) = mean(result.cpawl_all_NI(trange+1, 1, IDright), 'omitnan');
        pawguide_R_NI(2, i) = mean(result.cpawl_all_NI(trange+1, 2, IDright), 'omitnan');
        pawguide_R_NI(3, i) = mean(result.cpawl_all_NI(trange+1, 3, IDright), 'omitnan');
    end

    if ~(sum(result.IDleft_front_I-result.IDleft_top_I) == 0 && sum(result.IDright_front_I-result.IDright_top_I) == 0)
        errordlg(['Check ' filename{i} ' for the left and right IDs']);
        return;
    end
    IDleft = result.IDleft_front_I;
    IDright = result.IDright_front_I;
    IDall = [IDleft IDright];
    side_I(i, 1) = numel(IDleft);
    side_I(i, 2) = numel(IDright);
    
    if ~isempty(IDall)
        ctopI = squeeze(result.ctop_all_I(trange+1, :, IDall));
        cbottomI = squeeze(result.cbottom_all_I(trange+1, :, IDall));
        cpawlI = squeeze(result.cpawl_all_I(trange+1, :, IDall));
        cpawrI = squeeze(result.cpawr_all_I(trange+1, :, IDall));
        ceyelI = squeeze(result.ceyel_all_I(trange+1, 1, IDall));
        if numel(IDall) == 1
            ctopI = ctopI';
            cbottomI = cbottomI';
            cpawlI = cpawlI';
            cpawrI = cpawrI';
            ceyelI = ceyelI';
        end
    else
        ctopI = [];
        cbottomI = [];
        cpawlI = [];
        cpawrI = [];
        ceyelI = [];
    end
    
    if ~isempty(IDall)
        if numel(IDall) > thr
            ceyel_I(:, i, 1) = mean(squeeze(result.ceyel_all_I(trange+1, :, IDall)), 2, 'omitnan');
%             ceyer_I(:, i, 1) = mean(squeeze(result.ceyer_all_I(trange+1, :, IDall)), 2, 'omitnan');
            cnose_I(:, i, 1) = mean(squeeze(result.cnose_all_I(trange+1, :, IDall)), 2, 'omitnan');
            cpawl_I(:, i, 1) = mean(squeeze(result.cpawl_all_I(trange+1, :, IDall)), 2, 'omitnan');
            cpawr_I(:, i, 1) = mean(squeeze(result.cpawr_all_I(trange+1, :, IDall)), 2, 'omitnan');
            ctop_I(:, i, 1) = mean(squeeze(result.ctop_all_I(trange+1, :, IDall)), 2, 'omitnan');
            cbottom_I(:, i, 1) = mean(squeeze(result.cbottom_all_I(trange+1, :, IDall)), 2, 'omitnan');
        end
        if numel(IDleft) > thr
            ceyel_I(:, i, 2) = mean(squeeze(result.ceyel_all_I(trange+1, :, IDleft)), 2, 'omitnan');
            %             ceyer_I(:, i, 2) = mean(squeeze(result.ceyer_all_I(trange+1, :, IDleft)), 2, 'omitnan');
            cnose_I(:, i, 2) = mean(squeeze(result.cnose_all_I(trange+1, :, IDleft)), 2, 'omitnan');
            cpawl_I(:, i, 2) = mean(squeeze(result.cpawl_all_I(trange+1, :, IDleft)), 2, 'omitnan');
            cpawr_I(:, i, 2) = mean(squeeze(result.cpawr_all_I(trange+1, :, IDleft)), 2, 'omitnan');
            ctop_I(:, i, 2) = mean(squeeze(result.ctop_all_I(trange+1, :, IDleft)), 2, 'omitnan');
            cbottom_I(:, i, 2) = mean(squeeze(result.cbottom_all_I(trange+1, :, IDleft)), 2, 'omitnan');
            
%             pawgrasp_L_I(1, i) = mean(abs(result.cpawl_all_I(trange+1, 1, IDleft)-result.ceyel_all_I(trange+1, 1, IDleft)/2), 'omitnan');
            pawgrasp_L_I(1, i) = mean(result.cpawl_all_I(trange+1, 1, IDleft), 'omitnan');
            pawgrasp_L_I(2, i) = mean(result.cpawl_all_I(trange+1, 2, IDleft), 'omitnan');
            pawgrasp_L_I(3, i) = mean(result.cpawl_all_I(trange+1, 3, IDleft), 'omitnan');
%             pawguide_L_I(1, i) = mean(abs(result.cpawr_all_I(trange+1, 1, IDleft)-result.ceyel_all_I(trange+1, 1, IDleft)/2), 'omitnan');
            pawguide_L_I(1, i) = mean(result.cpawr_all_I(trange+1, 1, IDleft), 'omitnan');
            pawguide_L_I(2, i) = mean(result.cpawr_all_I(trange+1, 2, IDleft), 'omitnan');
            pawguide_L_I(3, i) = mean(result.cpawr_all_I(trange+1, 3, IDleft), 'omitnan');
        end
        if numel(IDright) > thr
            ceyel_I(:, i, 3) = mean(squeeze(result.ceyel_all_I(trange+1, :, IDright)), 2, 'omitnan');
%             ceyer_I(:, i, 3) = mean(squeeze(result.ceyer_all_I(trange+1, :, IDright)), 2, 'omitnan');
            cnose_I(:, i, 3) = mean(squeeze(result.cnose_all_I(trange+1, :, IDright)), 2, 'omitnan');
            cpawl_I(:, i, 3) = mean(squeeze(result.cpawl_all_I(trange+1, :, IDright)), 2, 'omitnan');
            cpawr_I(:, i, 3) = mean(squeeze(result.cpawr_all_I(trange+1, :, IDright)), 2, 'omitnan');
            ctop_I(:, i, 3) = mean(squeeze(result.ctop_all_I(trange+1, :, IDright)), 2, 'omitnan');
            cbottom_I(:, i, 3) = mean(squeeze(result.cbottom_all_I(trange+1, :, IDright)), 2, 'omitnan');
            
%             pawgrasp_R_I(1, i) = mean(abs(result.cpawr_all_I(trange+1, 1, IDright)-result.ceyel_all_I(trange+1, 1, IDright)/2), 'omitnan');
            pawgrasp_R_I(1, i) = mean(result.cpawr_all_I(trange+1, 1, IDright), 'omitnan');
            pawgrasp_R_I(2, i) = mean(result.cpawr_all_I(trange+1, 2, IDright), 'omitnan');
            pawgrasp_R_I(3, i) = mean(result.cpawr_all_I(trange+1, 3, IDright), 'omitnan');
%             pawguide_R_I(1, i) = mean(abs(result.cpawl_all_I(trange+1, 1, IDright)-result.ceyel_all_I(trange+1, 1, IDright)/2), 'omitnan');
            pawguide_R_I(1, i) = mean(result.cpawl_all_I(trange+1, 1, IDright), 'omitnan');
            pawguide_R_I(2, i) = mean(result.cpawl_all_I(trange+1, 2, IDright), 'omitnan');
            pawguide_R_I(3, i) = mean(result.cpawl_all_I(trange+1, 3, IDright), 'omitnan');
        end
    end
    
    xmidline(i) = mean([ceyelNI; ceyelI])/2;
    
    [orientationAH_NI, orientationAH_I] = compute_projected_orientation(ctopNI, cbottomNI, ctopI, cbottomI, pcheck, stp);
    if pcheck
        set(gcf, 'Name', 'Pasta Orientation');
    end
    [orientationPaw_NI, orientationPaw_I] = compute_projected_orientation(cpawlNI, cpawrNI, cpawlI, cpawrI, pcheck, stp);
    if pcheck
        set(gcf, 'Name', 'Two Paw Orientation');
    end
    orientation_all_NI.xy = [orientation_all_NI.xy orientationAH_NI.xy];
    orientation_all_NI.xz = [orientation_all_NI.xz orientationAH_NI.xz];
    orientation_all_NI.yz = [orientation_all_NI.yz orientationAH_NI.yz];
    orientation_all_NI.xy_avg(i) = mean(abs(orientationAH_NI.xy-90));
    orientation_all_NI.xz_avg(i) = mean(abs(orientationAH_NI.xz-90));
    orientation_all_NI.yz_avg(i) = mean(orientationAH_NI.yz);
    orientation_all_I.xy = [orientation_all_I.xy orientationAH_I.xy];
    orientation_all_I.xz = [orientation_all_I.xz orientationAH_I.xz];
    orientation_all_I.yz = [orientation_all_I.yz orientationAH_I.yz];
    if ~isempty(IDall)
        orientation_all_I.xy_avg(i) = mean(abs(orientationAH_I.xy-90));
        orientation_all_I.xz_avg(i) = mean(abs(orientationAH_I.xz-90));
        orientation_all_I.yz_avg(i) = mean(orientationAH_I.yz);
    end
    
%     orientationAH_NI.xy(orientationAH_NI.xy > 90) = 180-orientationAH_NI.xy(orientationAH_NI.xy > 90);
%     if ~isempty(orientationAH_I.xy)
%         orientationAH_I.xy(orientationAH_I.xy > 90) = 180-orientationAH_I.xy(orientationAH_I.xy > 90);
%     end
    
    orientationPaw_NI.xy(orientationPaw_NI.xy > 90) = orientationPaw_NI.xy(orientationPaw_NI.xy > 90)-180;
    if ~isempty(orientationPaw_I.xy)
        orientationPaw_I.xy(orientationPaw_I.xy > 90) = orientationPaw_I.xy(orientationPaw_I.xy > 90)-180;
    end
    orientationPaw_NI.xz(orientationPaw_NI.xz > 90) = orientationPaw_NI.xz(orientationPaw_NI.xz > 90)-180;
    if ~isempty(orientationPaw_I.xz)
        orientationPaw_I.xz(orientationPaw_I.xz > 90) = orientationPaw_I.xz(orientationPaw_I.xz > 90)-180;
    end
    
    figure;
    if ~isempty(cpawrI) && ~isempty(cpawlI)
        subplot(1, 3, 1);
        correlation_plot(cpawlNI(1, :), cpawrNI(1, :), side_NI(i, 1), cpawlI(1, :), cpawrI(1, :), side_I(i, 1), 'PawL X (mm)', 'PawR X (mm)', 0);
        subplot(1, 3, 2);
        correlation_plot(cpawlNI(2, :), cpawrNI(2, :), side_NI(i, 1), cpawlI(2, :), cpawrI(2, :), side_I(i, 1), 'PawL Y (mm)', 'PawR Y (mm)', 0);
        subplot(1, 3, 3);
        correlation_plot(cpawlNI(3, :), cpawrNI(3, :), side_NI(i, 1), cpawlI(3, :), cpawrI(3, :), side_I(i, 1), 'PawL Z (mm)', 'PawR Z (mm)', 0);
    else
        subplot(1, 3, 1);
        correlation_plot(cpawlNI(1, :), cpawrNI(1, :), side_NI(i, 1), [], [], [], 'PawL X (mm)', 'PawR X (mm)', 0);
        subplot(1, 3, 2);
        correlation_plot(cpawlNI(2, :), cpawrNI(2, :), side_NI(i, 1), [], [], [], 'PawL Y (mm)', 'PawR Y (mm)', 0);
        subplot(1, 3, 3);
        correlation_plot(cpawlNI(3, :), cpawrNI(3, :), side_NI(i, 1), [], [], [], 'PawL Z (mm)', 'PawR Z (mm)', 0);
    end
    
    figure;
    subplot(2, 4, 1);
    correlation_plot(orientationPaw_NI.xy, orientationAH_NI.xy, side_NI(i, 1), orientationPaw_I.xy, orientationAH_I.xy, side_I(i, 1), ['Paw XY Orientation (' char(176) ')'], ['AH XY Orientation (' char(176) ')'], 1);
    if ~isempty(cpawrI)
        subplot(2, 4, 2);
        result = correlation_plot(cpawrNI(1, :), orientationAH_NI.xy, side_NI(i, 1), cpawrI(1, :), orientationAH_I.xy, side_I(i, 1), 'PawR X (mm)', ['AH XY Orientation (' char(176) ')'], 1);
        corrX_pawr.rsquare(i, :) = result.rsquare;
        corrX_pawr.slope(i, :) = result.slope;
        corrX_pawr.rho(i, :) = result.rho;
        corrX_pawr.p(i, :) = result.p;
        subplot(2, 4, 3);
        result = correlation_plot(cpawrNI(2, :), orientationAH_NI.xy, side_NI(i, 1), cpawrI(2, :), orientationAH_I.xy, side_I(i, 1), 'PawR Y (mm)', ['AH XY Orientation (' char(176) ')'], 1);
        corrY_pawr.rsquare(i, :) = result.rsquare;
        corrY_pawr.slope(i, :) = result.slope;
        corrY_pawr.rho(i, :) = result.rho;
        corrY_pawr.p(i, :) = result.p;
        subplot(2, 4, 4);
        result = correlation_plot(cpawrNI(3, :), orientationAH_NI.xy, side_NI(i, 1), cpawrI(3, :), orientationAH_I.xy, side_I(i, 1), 'PawR Z (mm)', ['AH XY Orientation (' char(176) ')'], 1);
        corrZ_pawr.rsquare(i, :) = result.rsquare;
        corrZ_pawr.slope(i, :) = result.slope;
        corrZ_pawr.rho(i, :) = result.rho;
        corrZ_pawr.p(i, :) = result.p;
    else
        subplot(2, 4, 2);
        result = correlation_plot(cpawrNI(1, :), orientationAH_NI.xy, side_NI(i, 1), [], orientationAH_I.xy, side_I(i, 1), 'PawR X (mm)', ['AH XY Orientation (' char(176) ')'], 1);
        corrX_pawr.rsquare(i, :) = result.rsquare;
        corrX_pawr.slope(i, :) = result.slope;
        corrX_pawr.rho(i, :) = result.rho;
        corrX_pawr.p(i, :) = result.p;
        subplot(2, 4, 3);
        result = correlation_plot(cpawrNI(2, :), orientationAH_NI.xy, side_NI(i, 1), [], orientationAH_I.xy, side_I(i, 1), 'PawR Y (mm)', ['AH XY Orientation (' char(176) ')'], 1);
        corrY_pawr.rsquare(i, :) = result.rsquare;
        corrY_pawr.slope(i, :) = result.slope;
        corrY_pawr.rho(i, :) = result.rho;
        corrY_pawr.p(i, :) = result.p;
        subplot(2, 4, 4);
        result = correlation_plot(cpawrNI(3, :), orientationAH_NI.xy, side_NI(i, 1), [], orientationAH_I.xy, side_I(i, 1), 'PawR Z (mm)', ['AH XY Orientation (' char(176) ')'], 1);
        corrZ_pawr.rsquare(i, :) = result.rsquare;
        corrZ_pawr.slope(i, :) = result.slope;
        corrZ_pawr.rho(i, :) = result.rho;
        corrZ_pawr.p(i, :) = result.p;
    end
    subplot(2, 4, 5);
    correlation_plot(orientationPaw_NI.xz, orientationAH_NI.xy, side_NI(i, 1), orientationPaw_I.xz, orientationAH_I.xy, side_I(i, 1), ['Paw XZ Orientation (' char(176) ')'], ['AH XY Orientation (' char(176) ')'], 1);
    if ~isempty(cpawlI)
        subplot(2, 4, 6);
        result = correlation_plot(cpawlNI(1, :), orientationAH_NI.xy, side_NI(i, 1), cpawlI(1, :), orientationAH_I.xy, side_I(i, 1), 'PawL X (mm)', ['AH XY Orientation (' char(176) ')'], 1);
        corrX_pawl.rsquare(i, :) = result.rsquare;
        corrX_pawl.slope(i, :) = result.slope;
        corrX_pawl.rho(i, :) = result.rho;
        corrX_pawl.p(i, :) = result.p;
        subplot(2, 4, 7);
        result = correlation_plot(cpawlNI(2, :), orientationAH_NI.xy, side_NI(i, 1), cpawlI(2, :), orientationAH_I.xy, side_I(i, 1), 'PawL Y (mm)', ['AH XY Orientation (' char(176) ')'], 1);
        corrY_pawl.rsquare(i, :) = result.rsquare;
        corrY_pawl.slope(i, :) = result.slope;
        corrY_pawl.rho(i, :) = result.rho;
        corrY_pawl.p(i, :) = result.p;
        subplot(2, 4, 8);
        result = correlation_plot(cpawlNI(3, :), orientationAH_NI.xy, side_NI(i, 1), cpawlI(3, :), orientationAH_I.xy, side_I(i, 1), 'PawL Z (mm)', ['AH XY Orientation (' char(176) ')'], 1);
        corrZ_pawl.rsquare(i, :) = result.rsquare;
        corrZ_pawl.slope(i, :) = result.slope;
        corrZ_pawl.rho(i, :) = result.rho;
        corrZ_pawl.p(i, :) = result.p;
    else
        subplot(2, 4, 6);
        result = correlation_plot(cpawlNI(1, :), orientationAH_NI.xy, side_NI(i, 1), [], orientationAH_I.xy, side_I(i, 1), 'PawL X (mm)', ['AH XY Orientation (' char(176) ')'], 1);
        corrX_pawl.rsquare(i, :) = result.rsquare;
        corrX_pawl.slope(i, :) = result.slope;
        corrX_pawl.rho(i, :) = result.rho;
        corrX_pawl.p(i, :) = result.p;
        subplot(2, 4, 7);
        result = correlation_plot(cpawlNI(2, :), orientationAH_NI.xy, side_NI(i, 1), [], orientationAH_I.xy, side_I(i, 1), 'PawL Y (mm)', ['AH XY Orientation (' char(176) ')'], 1);
        corrY_pawl.rsquare(i, :) = result.rsquare;
        corrY_pawl.slope(i, :) = result.slope;
        corrY_pawl.rho(i, :) = result.rho;
        corrY_pawl.p(i, :) = result.p;
        subplot(2, 4, 8);
        result = correlation_plot(cpawlNI(3, :), orientationAH_NI.xy, side_NI(i, 1), [], orientationAH_I.xy, side_I(i, 1), 'PawL Z (mm)', ['AH XY Orientation (' char(176) ')'], 1);
        corrZ_pawl.rsquare(i, :) = result.rsquare;
        corrZ_pawl.slope(i, :) = result.slope;
        corrZ_pawl.rho(i, :) = result.rho;
        corrZ_pawl.p(i, :) = result.p;
    end
end
figure;
subplot(3, 1, 1);
% draw_orientations(abs(orientation_all_NI.xz-90), abs(orientation_all_I.xz-90), edges, 'Orientation XZ');
draw_orientations(orientation_all_NI.xz, orientation_all_I.xz, edges, 'Orientation XZ');
subplot(3, 1, 2);
% draw_orientations(abs(orientation_all_NI.xy-90), abs(orientation_all_I.xy-90), edges, 'Orientation XY');
draw_orientations(orientation_all_NI.xy, orientation_all_I.xy, edges, 'Orientation XY');
subplot(3, 1, 3);
draw_orientations(orientation_all_NI.yz, orientation_all_I.yz, edges, 'Orientation YZ');

figure;
disp('Orientation XY');
subplot(1, 3, 1);
paired_plot(orientation_all_NI.xy_avg, orientation_all_I.xy_avg, ['AH XY Orientation (' char(176) ')'], filename, 'bar');
disp('Orientation XZ');
subplot(1, 3, 2);
paired_plot(orientation_all_NI.xz_avg, orientation_all_I.xz_avg, ['AH XZ Orientation (' char(176) ')'], filename, 'bar');
disp('Orientation YZ');
subplot(1, 3, 3);
paired_plot(orientation_all_NI.yz_avg, orientation_all_I.yz_avg, ['AH YZ Orientation (' char(176) ')'], filename, 'bar');

figure;
disp('Grasp paw X to midline');
subplot(2, 3, 1);
paired_plot(abs([pawgrasp_L_NI(1, :) pawgrasp_R_NI(1, :)]-[xmidline xmidline]), abs([pawgrasp_L_I(1, :) pawgrasp_R_I(1, :)]-[xmidline xmidline]), 'Grasp paw X to midline (mm)', [filename filename], 'dot');
disp('Grasp paw Y');
subplot(2, 3, 2);
paired_plot([pawgrasp_L_NI(2, :) pawgrasp_R_NI(2, :)], [pawgrasp_L_I(2, :) pawgrasp_R_I(2, :)], 'Grasp paw Y (mm)', [filename filename], 'dot');
disp('Grasp paw Z');
subplot(2, 3, 3);
paired_plot([pawgrasp_L_NI(3, :) pawgrasp_R_NI(3, :)], [pawgrasp_L_I(3, :) pawgrasp_R_I(3, :)], 'Grasp paw Z (mm)', [filename filename], 'dot');
disp('Guide paw X to midline');
subplot(2, 3, 4);
paired_plot(abs([pawguide_L_NI(1, :) pawguide_R_NI(1, :)]-[xmidline xmidline]), abs([pawguide_L_I(1, :) pawguide_R_I(1, :)]-[xmidline xmidline]), 'Guide paw X to midline (mm)', [filename filename], 'dot');
disp('Guide paw Y');
subplot(2, 3, 5);
paired_plot([pawguide_L_NI(2, :) pawguide_R_NI(2, :)], [pawguide_L_I(2, :) pawguide_R_I(2, :)], 'Guide paw Y (mm)', [filename filename], 'dot');
disp('Guide paw Z');
subplot(2, 3, 6);
paired_plot([pawguide_L_NI(3, :) pawguide_R_NI(3, :)], [pawguide_L_I(3, :) pawguide_R_I(3, :)], 'Guide paw Z (mm)', [filename filename], 'dot');

figure;
subplot(3, 1, 1);
disp('Grasp paw X to midline');
paired_plot(abs([pawgrasp_L_NI(1, :) pawgrasp_R_NI(1, :)]-[xmidline xmidline]), abs([pawgrasp_L_I(1, :) pawgrasp_R_I(1, :)]-[xmidline xmidline]), 'Grasp paw X to midline (mm)', [filename filename], 'dot', [1 2]);
disp('Guide paw X to midline');
paired_plot(abs([pawguide_L_NI(1, :) pawguide_R_NI(1, :)]-[xmidline xmidline]), abs([pawguide_L_I(1, :) pawguide_R_I(1, :)]-[xmidline xmidline]), 'Guide paw X to midline (mm)', [filename filename], 'dot', [3 4]);
xlim([0.4 4.6]);
subplot(3, 1, 2);
disp('Grasp paw Y');
paired_plot([pawgrasp_L_NI(2, :) pawgrasp_R_NI(2, :)], [pawgrasp_L_I(2, :) pawgrasp_R_I(2, :)], 'Grasp paw Y (mm)', [filename filename], 'dot', [1 2]);
disp('Guide paw Y');
paired_plot([pawguide_L_NI(2, :) pawguide_R_NI(2, :)], [pawguide_L_I(2, :) pawguide_R_I(2, :)], 'Guide paw Y (mm)', [filename filename], 'dot', [3 4]);
xlim([0.4 4.6]);
subplot(3, 1, 3);
disp('Grasp paw Z');
paired_plot([pawgrasp_L_NI(3, :) pawgrasp_R_NI(3, :)], [pawgrasp_L_I(3, :) pawgrasp_R_I(3, :)], 'Grasp paw Z (mm)', [filename filename], 'dot', [1 2]);
disp('Guide paw Z');
paired_plot([pawguide_L_NI(3, :) pawguide_R_NI(3, :)], [pawguide_L_I(3, :) pawguide_R_I(3, :)], 'Guide paw Z (mm)', [filename filename], 'dot', [3 4]);
xlim([0.4 4.6]);

title_text = {'All', 'Left', 'Right'};
for i = 1:3
%     figure('Name', title_text{i});
%     subplot(2, 3, 1);
%     paired_plot(cpawl_NI(1, :, i), cpawl_I(1, :, i), 'X (mm)', filename, 'dot');
%     subplot(2, 3, 2);
%     paired_plot(cpawl_NI(2, :, i), cpawl_I(2, :, i), 'Y (mm)', filename, 'dot');
%     subplot(2, 3, 3);
%     paired_plot(cpawl_NI(3, :, i), cpawl_I(3, :, i), 'Z (mm)', filename, 'dot');
%     subplot(2, 3, 4);
%     paired_plot(cpawr_NI(1, :, i), cpawr_I(1, :, i), 'X (mm)', filename, 'dot');
%     subplot(2, 3, 5);
%     paired_plot(cpawr_NI(2, :, i), cpawr_I(2, :, i), 'Y (mm)', filename, 'dot');
%     subplot(2, 3, 6);
%     paired_plot(cpawr_NI(3, :, i), cpawr_I(3, :, i), 'Z (mm)', filename, 'dot');
%     
%     draw_individual_posture(cnose_NI(:, :, i), ceyel_NI(:, :, i), ceyer_NI(:, :, i), cpawl_NI(:, :, i), cpawr_NI(:, :, i), ctop_NI(:, :, i), cbottom_NI(:, :, i), bitelocation_NI(i, :, :),...
%         cnose_I(:, :, i), ceyel_I(:, :, i), ceyer_I(:, :, i), cpawl_I(:, :, i), cpawr_I(:, :, i), ctop_I(:, :, i), cbottom_I(:, :, i), bitelocation_I(i, :, :), filename);
%     set(gcf, 'Name', title_text{i});
    
    ref = ceyer_NI(:, :, i);
    cnose_NI(:, :, i) = cnose_NI(:, :, i)-ref;
    ceyel_NI(:, :, i) = ceyel_NI(:, :, i)-ref;
    ceyer_NI(:, :, i) = ceyer_NI(:, :, i)-ref;
    cpawl_NI(:, :, i) = cpawl_NI(:, :, i)-ref;
    cpawr_NI(:, :, i) = cpawr_NI(:, :, i)-ref;
    ctop_NI(:, :, i) = ctop_NI(:, :, i)-ref;
    cbottom_NI(:, :, i) = cbottom_NI(:, :, i)-ref;
    
    ref = ceyer_I(:, :, i);
    cnose_I(:, :, i) = cnose_I(:, :, i)-ref;
    ceyel_I(:, :, i) = ceyel_I(:, :, i)-ref;
    ceyer_I(:, :, i) = ceyer_I(:, :, i)-ref;
    cpawl_I(:, :, i) = cpawl_I(:, :, i)-ref;
    cpawr_I(:, :, i) = cpawr_I(:, :, i)-ref;
    ctop_I(:, :, i) = ctop_I(:, :, i)-ref;
    cbottom_I(:, :, i) = cbottom_I(:, :, i)-ref;
    
    figure('Name', title_text{i});
    disp('Left Paw X');
    subplot(2, 3, 1);
    paired_plot(cpawl_NI(1, :, i), cpawl_I(1, :, i), 'PawL X (mm)', filename, 'dot');
    disp('Left Paw Y');
    subplot(2, 3, 2);
    paired_plot(cpawl_NI(2, :, i), cpawl_I(2, :, i), 'PawL Y (mm)', filename, 'dot');
    disp('Left Paw Z');
    subplot(2, 3, 3);
    paired_plot(cpawl_NI(3, :, i), cpawl_I(3, :, i), 'PawL Z (mm)', filename, 'dot');
    disp('Right Paw X');
    subplot(2, 3, 4);
    paired_plot(cpawr_NI(1, :, i), cpawr_I(1, :, i), 'PawR X (mm)', filename, 'dot');
    disp('Right Paw Y');
    subplot(2, 3, 5);
    paired_plot(cpawr_NI(2, :, i), cpawr_I(2, :, i), 'PawR Y (mm)', filename, 'dot');
    disp('Right Paw Z');
    subplot(2, 3, 6);
    paired_plot(cpawr_NI(3, :, i), cpawr_I(3, :, i), 'PawR Z (mm)', filename, 'dot');
    
    draw_individual_posture(cnose_NI(:, :, i), ceyel_NI(:, :, i), ceyer_NI(:, :, i), cpawl_NI(:, :, i), cpawr_NI(:, :, i), ctop_NI(:, :, i), cbottom_NI(:, :, i),...
        cnose_I(:, :, i), ceyel_I(:, :, i), ceyer_I(:, :, i), cpawl_I(:, :, i), cpawr_I(:, :, i), ctop_I(:, :, i), cbottom_I(:, :, i),filename);
    set(gcf, 'Name', title_text{i});
end

cd(curpwd);

function draw_individual_posture(cnose_all_NI, ceyel_all_NI, ceyer_all_NI, cpawl_all_NI, cpawr_all_NI, ctop_all_NI, cbottom_all_NI,...
    cnose_all_I, ceyel_all_I, ceyer_all_I, cpawl_all_I, cpawr_all_I, ctop_all_I, cbottom_all_I, filename)
markers = {'o', 'x', '+', '*', 's', 'd', 'p', 'h', 'v', '^', '<', '>', '.'};
mksize = 6;
figure;
subplot(1, 2, 1);
hold on;
for j = 1:size(ceyel_all_NI, 2)
    plot([ceyer_all_NI(1, j) ceyel_all_NI(1, j)], [ceyer_all_NI(3, j) ceyel_all_NI(3, j)], '-k', 'ButtonDownFcn', @displayfilename, 'UserData', filename{j});
    plot([cnose_all_NI(1, j) ceyel_all_NI(1, j)], [cnose_all_NI(3, j) ceyel_all_NI(3, j)], '-k', 'ButtonDownFcn', @displayfilename, 'UserData', filename{j});
    plot([cnose_all_NI(1, j) ceyer_all_NI(1, j)], [cnose_all_NI(3, j) ceyer_all_NI(3, j)], '-k', 'ButtonDownFcn', @displayfilename, 'UserData', filename{j});
    plot(ceyel_all_NI(1, j), ceyel_all_NI(3, j), markers{j}, 'Color', [0 0 0], 'MarkerSize', mksize);
    plot(cnose_all_NI(1, j), cnose_all_NI(3, j), markers{j}, 'Color', [0.5 0 0], 'MarkerSize', mksize);
    
    plot([cnose_all_NI(1, j) cpawl_all_NI(1, j)], [cnose_all_NI(3, j) cpawl_all_NI(3, j)], '-k', 'ButtonDownFcn', @displayfilename, 'UserData', filename{j});
    plot([cnose_all_NI(1, j) cpawr_all_NI(1, j)], [cnose_all_NI(3, j) cpawr_all_NI(3, j)], '-k', 'ButtonDownFcn', @displayfilename, 'UserData', filename{j});
    plot([cpawl_all_NI(1, j) cpawr_all_NI(1, j)], [cpawl_all_NI(3, j) cpawr_all_NI(3, j)], '-k', 'ButtonDownFcn', @displayfilename, 'UserData', filename{j});
    
    plot(cpawl_all_NI(1, j), cpawl_all_NI(3, j), markers{j}, 'Color', [0.00,0.45,0.74], 'MarkerSize', mksize);
    plot(cpawr_all_NI(1, j), cpawr_all_NI(3, j), markers{j}, 'Color', [0.47,0.67,0.19], 'MarkerSize', mksize);
    
    plot([ctop_all_NI(1, j) cbottom_all_NI(1, j)], [ctop_all_NI(3, j) cbottom_all_NI(3, j)], '-',...
        'Marker', markers{j}, 'Color', [0.93,0.69,0.13], 'MarkerSize', mksize, 'ButtonDownFcn', @displayfilename, 'UserData', filename{j});
end

for j = 1:size(ceyel_all_I, 2)
    plot([ceyer_all_I(1, j) ceyel_all_I(1, j)], [ceyer_all_I(3, j) ceyel_all_I(3, j)], '-', 'Color', [0.8 0.8 0.8], 'ButtonDownFcn', @displayfilename, 'UserData', filename{j});
    plot([cnose_all_I(1, j) ceyel_all_I(1, j)], [cnose_all_I(3, j) ceyel_all_I(3, j)], '-', 'Color', [0.8 0.8 0.8], 'ButtonDownFcn', @displayfilename, 'UserData', filename{j});
    plot([cnose_all_I(1, j) ceyer_all_I(1, j)], [cnose_all_I(3, j) ceyer_all_I(3, j)], '-', 'Color', [0.8 0.8 0.8], 'ButtonDownFcn', @displayfilename, 'UserData', filename{j});
    plot(ceyel_all_I(1, j), ceyel_all_I(3, j), markers{j}, 'Color', [0.8 0.8 0.8], 'MarkerSize', mksize);
    plot(cnose_all_I(1, j), cnose_all_I(3, j), markers{j}, 'Color', [1 0 0], 'MarkerSize', mksize);
    
    plot([cnose_all_I(1, j) cpawl_all_I(1, j)], [cnose_all_I(3, j) cpawl_all_I(3, j)], '-', 'Color', [0.8 0.8 0.8], 'ButtonDownFcn', @displayfilename, 'UserData', filename{j});
    plot([cnose_all_I(1, j) cpawr_all_I(1, j)], [cnose_all_I(3, j) cpawr_all_I(3, j)], '-', 'Color', [0.8 0.8 0.8], 'ButtonDownFcn', @displayfilename, 'UserData', filename{j});
    plot([cpawl_all_I(1, j) cpawr_all_I(1, j)], [cpawl_all_I(3, j) cpawr_all_I(3, j)], '-',...
        'Color', [0.8 0.8 0.8], 'ButtonDownFcn', @displayfilename, 'UserData', filename{j});
    
    plot(cpawl_all_I(1, j), cpawl_all_I(3, j), markers{j}, 'Color', [0.30,0.75,0.93], 'MarkerSize', mksize);
    plot(cpawr_all_I(1, j), cpawr_all_I(3, j), markers{j}, 'Color', [0 1 0], 'MarkerSize', mksize);
    
    plot([ctop_all_I(1, j) cbottom_all_I(1, j)], [ctop_all_I(3, j) cbottom_all_I(3, j)], '-',...
        'Marker', markers{j}, 'Color', [0.92,0.92,0.45], 'MarkerSize', mksize, 'ButtonDownFcn', @displayfilename, 'UserData', filename{j});
end

axis equal;
xlabel('X (mm)');
ylabel('Z (mm)');
set(gca, 'TickDir', 'out', 'FontSize', 12, 'ButtonDownFcn', @extract_figure);

subplot(1, 2, 2);
hold on;
for j = 1:size(ceyel_all_NI, 2)
    plot([ceyer_all_NI(1, j) ceyel_all_NI(1, j)], [ceyer_all_NI(2, j) ceyel_all_NI(2, j)], '-k', 'ButtonDownFcn', @displayfilename, 'UserData', filename{j});
    plot([cnose_all_NI(1, j) ceyel_all_NI(1, j)], [cnose_all_NI(2, j) ceyel_all_NI(2, j)], '-k', 'ButtonDownFcn', @displayfilename, 'UserData', filename{j});
    plot([cnose_all_NI(1, j) ceyer_all_NI(1, j)], [cnose_all_NI(2, j) ceyer_all_NI(2, j)], '-k', 'ButtonDownFcn', @displayfilename, 'UserData', filename{j});
    plot(ceyel_all_NI(1, j), ceyel_all_NI(2, j), markers{j}, 'Color', [0 0 0], 'MarkerSize', mksize);
    plot(cnose_all_NI(1, j), cnose_all_NI(2, j), markers{j}, 'Color', [0.5 0 0], 'MarkerSize', mksize);
    
    plot([cnose_all_NI(1, j) cpawl_all_NI(1, j)], [cnose_all_NI(2, j) cpawl_all_NI(2, j)], '-k', 'ButtonDownFcn', @displayfilename, 'UserData', filename{j});
    plot([cnose_all_NI(1, j) cpawr_all_NI(1, j)], [cnose_all_NI(2, j) cpawr_all_NI(2, j)], '-k', 'ButtonDownFcn', @displayfilename, 'UserData', filename{j});
    plot([cpawl_all_NI(1, j) cpawr_all_NI(1, j)], [cpawl_all_NI(2, j) cpawr_all_NI(2, j)], '-k', 'ButtonDownFcn', @displayfilename, 'UserData', filename{j});
    
    plot(cpawl_all_NI(1, j), cpawl_all_NI(2, j), markers{j}, 'Color', [0.00,0.45,0.74], 'MarkerSize', mksize);
    plot(cpawr_all_NI(1, j), cpawr_all_NI(2, j), markers{j}, 'Color', [0.47,0.67,0.19], 'MarkerSize', mksize);
    
    plot([ctop_all_NI(1, j) cbottom_all_NI(1, j)], [ctop_all_NI(2, j) cbottom_all_NI(2, j)], '-',...
        'Marker', markers{j}, 'Color', [0.93,0.69,0.13], 'MarkerSize', mksize, 'ButtonDownFcn', @displayfilename, 'UserData', filename{j});
end

for j = 1:size(ceyel_all_I, 2)
    plot([ceyer_all_I(1, j) ceyel_all_I(1, j)], [ceyer_all_I(2, j) ceyel_all_I(2, j)], '-', 'Color', [0.8 0.8 0.8], 'ButtonDownFcn', @displayfilename, 'UserData', filename{j});
    plot([cnose_all_I(1, j) ceyel_all_I(1, j)], [cnose_all_I(2, j) ceyel_all_I(2, j)], '-', 'Color', [0.8 0.8 0.8], 'ButtonDownFcn', @displayfilename, 'UserData', filename{j});
    plot([cnose_all_I(1, j) ceyer_all_I(1, j)], [cnose_all_I(2, j) ceyer_all_I(2, j)], '-', 'Color', [0.8 0.8 0.8], 'ButtonDownFcn', @displayfilename, 'UserData', filename{j});
    plot(ceyel_all_I(1, j), ceyel_all_I(2, j), markers{j}, 'Color', [0.8 0.8 0.8], 'MarkerSize', mksize);
    plot(cnose_all_I(1, j), cnose_all_I(2, j), markers{j}, 'Color', [1 0 0], 'MarkerSize', mksize);
    
    plot([cnose_all_I(1, j) cpawl_all_I(1, j)], [cnose_all_I(2, j) cpawl_all_I(2, j)], '-', 'Color', [0.8 0.8 0.8], 'ButtonDownFcn', @displayfilename, 'UserData', filename{j});
    plot([cnose_all_I(1, j) cpawr_all_I(1, j)], [cnose_all_I(2, j) cpawr_all_I(2, j)], '-', 'Color', [0.8 0.8 0.8], 'ButtonDownFcn', @displayfilename, 'UserData', filename{j});
    plot([cpawl_all_I(1, j) cpawr_all_I(1, j)], [cpawl_all_I(2, j) cpawr_all_I(2, j)], '-',...
        'Color', [0.8 0.8 0.8], 'ButtonDownFcn', @displayfilename, 'UserData', filename{j});
    
    plot(cpawl_all_I(1, j), cpawl_all_I(2, j), markers{j}, 'Color', [0.30,0.75,0.93], 'MarkerSize', mksize);
    plot(cpawr_all_I(1, j), cpawr_all_I(2, j), markers{j}, 'Color', [0 1 0], 'MarkerSize', mksize);
    
    plot([ctop_all_I(1, j) cbottom_all_I(1, j)], [ctop_all_I(2, j) cbottom_all_I(2, j)], '-',...
        'Marker', markers{j}, 'Color', [0.92,0.92,0.45], 'MarkerSize', mksize, 'ButtonDownFcn', @displayfilename, 'UserData', filename{j});
end
axis equal;
xlabel('X (mm)');
ylabel('Y (mm)');
set(gca, 'TickDir', 'out', 'FontSize', 12, 'ButtonDownFcn', @extract_figure);

    function displayfilename(src, eventdata)
        htext = text(eventdata.IntersectionPoint(1), eventdata.IntersectionPoint(2), src.UserData, 'FontSize', 12, 'Color', 'r');
        pause(2);
        try
            delete(htext)
            clear htext;
        end
    end
end