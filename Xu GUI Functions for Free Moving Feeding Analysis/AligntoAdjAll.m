function AligntoAdjAll(app, Exp_Path, FrameRate, nmedian)
persistent path
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

ceyel_all_NI = [];
ceyer_all_NI = [];
cnose_all_NI = [];
cpawl_all_NI = [];
cpawr_all_NI = [];
ctop_all_NI = [];
cbottom_all_NI = [];
ccenter_all_NI = [];
ceyel_all_I = [];
ceyer_all_I = [];
cnose_all_I = [];
cpawl_all_I = [];
cpawr_all_I = [];
ctop_all_I = [];
cbottom_all_I = [];
ccenter_all_I = [];

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
FrameRate = round(FrameRate);
trange = app.bitesEditField.Value; % time range for bite alignment
trange = trange*FrameRate;
analyze_method = 1; % 1 for nose and two eyes plane transformation; 2 for simple nose and two paws trajectory; 3 for axis rotation based on two ankles

nadj_NI = zeros(1, ninhibition);
nadj_I = zeros(1, ninhibition);
adjinfo_NI = cell(0);
adjinfo_I = cell(0);

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
            LabelledEvents = [];
            temp = load([Exp_Paths{i} '\LabelledEvents' num2str(j) '.mat']);
            LabelledEvents = temp.LabelledEvents;
            PawRAdjustmentStart = LabelledEvents.PawRAdjustmentStart;
            PawLAdjustmentStart = LabelledEvents.PawLAdjustmentStart;
            PawRAdjustmentEnd = LabelledEvents.PawRAdjustmentEnd;
            PawLAdjustmentEnd = LabelledEvents.PawLAdjustmentEnd;
            adjustment_start = get_adjustment_start(PawRAdjustmentStart, PawRAdjustmentEnd, PawLAdjustmentStart, PawLAdjustmentEnd);
            if isempty(adjustment_start)
                continue;
            end
            
            if analyze_method == 2
                [~, ~, ~, ~, ~, xnose, ynose, znose, speednose, accelerationnose, laserstart, laserstop] = trajectory_postprocessing(3, Exp_Paths{i}, j, label_table, nmedian, FrameRate);
                [~, ~, ~, ~, ~, xpawl, ypawl, zpawl, speedpawl, accelerationpawl, ~, ~] = trajectory_postprocessing(10, Exp_Paths{i}, j, label_table, nmedian, FrameRate);
                [~, ~, ~, ~, ~, xpawr, ypawr, zpawr, speedpawr, accelerationpawr, ~, ~] = trajectory_postprocessing(16, Exp_Paths{i}, j, label_table, nmedian, FrameRate);
                
                t = (1:size(xnose, 1))'/FrameRate;
                
                if app.NoLightButton.Value
                    laserstart = [];
                    laserstop = [];
                else
                    laser_timestamps = Bite_events(j).laser_timestamps;
                end
                
                if ~isempty(laserstart)
                    if ~isempty(adjustment_start)
                        for k = 1:ninhibition
                            for m = 1:numel(adjustment_start)
                                timestamp = adjustment_start(m);
                                if timestamp >= laser_timestamps(2*k-1) && timestamp <= laser_timestamps(2*k)
                                    [~, id] = min(abs(t-timestamp));
                                    if id-trange >= 1
                                        nadj_I(k) = nadj_I(k)+1;
                                        znoseadj_I(:, nadj_I(k), k) = znose(id-trange:id+trange, 1);
                                        zspeednoseadj_I(:, nadj_I(k), k) = abs(speednose(id-trange:id+trange, 3));
                                        xyzspeednoseadj_I(:, nadj_I(k), k) = speednose(id-trange:id+trange, 7);
                                        zaccelerationnoseadj_I(:, nadj_I(k), k) = abs(accelerationnose(id-trange:id+trange, 3));
                                        xyzaccelerationnoseadj_I(:, nadj_I(k), k) = abs(accelerationnose(id-trange:id+trange, 7));
                                        
                                        zpawladj_I(:, nadj_I(k), k) = zpawl(id-trange:id+trange, 1);
                                        zspeedpawladj_I(:, nadj_I(k), k) = abs(speedpawl(id-trange:id+trange, 3));
                                        xyzspeedpawladj_I(:, nadj_I(k), k) = speedpawl(id-trange:id+trange, 7);
                                        zaccelerationpawladj_I(:, nadj_I(k), k) = abs(accelerationpawl(id-trange:id+trange, 3));
                                        xyzaccelerationpawladj_I(:, nadj_I(k), k) = abs(accelerationpawl(id-trange:id+trange, 7));
                                        
                                        zpawradj_I(:, nadj_I(k), k) = zpawr(id-trange:id+trange, 1);
                                        zspeedpawradj_I(:, nadj_I(k), k) = abs(speedpawr(id-trange:id+trange, 3));
                                        xyzspeedpawradj_I(:, nadj_I(k), k) = speedpawr(id-trange:id+trange, 7);
                                        zaccelerationpawradj_I(:, nadj_I(k), k) = abs(accelerationpawr(id-trange:id+trange, 3));
                                        xyzaccelerationpawradj_I(:, nadj_I(k), k) = abs(accelerationpawr(id-trange:id+trange, 7));
                                    end
                                end
                            end
                        end
                    end
                else
                    if ~isempty(adjustment_start)
                        for k = 1:ninhibition
                            for m = 1:numel(adjustment_start)
                                timestamp = adjustment_start(m);
                                if timestamp >= tinhibition(2*k-1) && timestamp <= tinhibition(2*k)
                                    [~, id] = min(abs(t-timestamp));
                                    if id-trange >= 1
                                        nadj_NI(k) = nadj_NI(k)+1;
                                        znoseadj_NI(:, nadj_NI(k), k) = znose(id-trange:id+trange, 1);
                                        zspeednoseadj_NI(:, nadj_NI(k), k) = abs(speednose(id-trange:id+trange, 3));
                                        xyzspeednoseadj_NI(:, nadj_NI(k), k) = speednose(id-trange:id+trange, 7);
                                        zaccelerationnoseadj_NI(:, nadj_NI(k), k) = abs(accelerationnose(id-trange:id+trange, 3));
                                        xyzaccelerationnoseadj_NI(:, nadj_NI(k), k) = abs(accelerationnose(id-trange:id+trange, 7));
                                        
                                        zpawladj_NI(:, nadj_NI(k), k) = zpawl(id-trange:id+trange, 1);
                                        zspeedpawladj_NI(:, nadj_NI(k), k) = abs(speedpawl(id-trange:id+trange, 3));
                                        xyzspeedpawladj_NI(:, nadj_NI(k), k) = speedpawl(id-trange:id+trange, 7);
                                        zaccelerationpawladj_NI(:, nadj_NI(k), k) = abs(accelerationpawl(id-trange:id+trange, 3));
                                        xyzaccelerationpawladj_NI(:, nadj_NI(k), k) = abs(accelerationpawl(id-trange:id+trange, 7));
                                        
                                        zpawradj_NI(:, nadj_NI(k), k) = zpawr(id-trange:id+trange, 1);
                                        zspeedpawradj_NI(:, nadj_NI(k), k) = abs(speedpawr(id-trange:id+trange, 3));
                                        xyzspeedpawradj_NI(:, nadj_NI(k), k) = speedpawr(id-trange:id+trange, 7);
                                        zaccelerationpawradj_NI(:, nadj_NI(k), k) = abs(accelerationpawr(id-trange:id+trange, 3));
                                        xyzaccelerationpawradj_NI(:, nadj_NI(k), k) = abs(accelerationpawr(id-trange:id+trange, 7));
                                    end
                                end
                            end
                        end
                    end
                end
            elseif analyze_method == 1 || analyze_method == 3
                [~, ~, ~, ~, ~, xnose, ynose, znose, ~, ~, laserstart, laserstop] = trajectory_postprocessing(3, Exp_Paths{i}, j, label_table, nmedian, FrameRate);
                [~, ~, ~, ~, ~, xpawl, ypawl, zpawl, ~, ~, ~, ~] = trajectory_postprocessing(10, Exp_Paths{i}, j, label_table, nmedian, FrameRate);
                [~, ~, ~, ~, ~, xpawr, ypawr, zpawr, ~, ~, ~, ~] = trajectory_postprocessing(16, Exp_Paths{i}, j, label_table, nmedian, FrameRate);
                [y_topPG1, y_topPG3, z_topPG1, ~, z_topPG3, xtop, ytop, ztop, ~, ~, ~, ~] = trajectory_postprocessing(21, Exp_Paths{i}, j, label_table, nmedian, FrameRate);
                [y_bottomPG1, y_bottomPG3, z_bottomPG1, ~, z_bottomPG3, xbottom, ybottom, zbottom, ~, ~, ~, ~] = trajectory_postprocessing(22, Exp_Paths{i}, j, label_table, nmedian, FrameRate);
                [y_centerPG1, y_centerPG3, z_centerPG1, ~, z_centerPG3, xcenter, ycenter, zcenter, ~, ~, ~, ~] = trajectory_postprocessing(33, Exp_Paths{i}, j, label_table, nmedian, FrameRate);
                
                if analyze_method == 1
                    [~, ~, ~, ~, ~, xeyel, yeyel, zeyel, ~, ~, ~, ~] = trajectory_postprocessing(1, Exp_Paths{i}, j, label_table, nmedian, FrameRate);
                    [~, ~, ~, ~, ~, xeyer, yeyer, zeyer, ~, ~, ~, ~] = trajectory_postprocessing(2, Exp_Paths{i}, j, label_table, nmedian, FrameRate);
%                     corigin = [xnose(:, 1) ynose(:, 1) znose(:, 1)];
%                     % translation, new coordinate origin is nose
                    corigin = [xeyer(:, 1) yeyer(:, 1) zeyer(:, 1)];
                    % translation, new coordinate origin is right eye
                    ceyelT = [xeyel(:, 1)-corigin(:, 1) yeyel(:, 1)-corigin(:, 2) zeyel(:, 1)-corigin(:, 3)];
%                     ceyerT = [xeyer(:, 1)-corigin(:, 1) yeyer(:, 1)-corigin(:, 2) zeyer(:, 1)-corigin(:, 3)];
                    cnoseT = [xnose(:, 1)-corigin(:, 1) ynose(:, 1)-corigin(:, 2) znose(:, 1)-corigin(:, 3)];
                    cpawlT = [xpawl(:, 1)-corigin(:, 1) ypawl(:, 1)-corigin(:, 2) zpawl(:, 1)-corigin(:, 3)];
                    cpawrT = [xpawr(:, 1)-corigin(:, 1) ypawr(:, 1)-corigin(:, 2) zpawr(:, 1)-corigin(:, 3)];
                    ctopT = [xtop(:, 1)-corigin(:, 1) ytop(:, 1)-corigin(:, 2) ztop(:, 1)-corigin(:, 3)];
                    cbottomT = [xbottom(:, 1)-corigin(:, 1) ybottom(:, 1)-corigin(:, 2) zbottom(:, 1)-corigin(:, 3)];
                    ccenterT = [xcenter(:, 1)-corigin(:, 1) ycenter(:, 1)-corigin(:, 2) zcenter(:, 1)-corigin(:, 3)];
                    
                    y_topPG1T = y_topPG1-corigin(:, 2);
                    y_topPG3T = y_topPG3-corigin(:, 2);
                    z_topPG1T = z_topPG1-corigin(:, 3);
                    z_topPG3T = z_topPG3-corigin(:, 3);
                    y_bottomPG1T = y_bottomPG1-corigin(:, 2);
                    y_bottomPG3T = y_bottomPG3-corigin(:, 2);
                    z_bottomPG1T = z_bottomPG1-corigin(:, 3);
                    z_bottomPG3T = z_bottomPG3-corigin(:, 3);
                    y_centerPG1T = y_centerPG1-corigin(:, 2);
                    y_centerPG3T = y_centerPG3-corigin(:, 2);
                    z_centerPG1T = z_centerPG1-corigin(:, 3);
                    z_centerPG3T = z_centerPG3-corigin(:, 3);
                    y_pastaT = [y_topPG1T y_bottomPG1T y_centerPG1T y_topPG3T y_bottomPG3T y_centerPG3T];
                    z_pastaT = [z_topPG1T z_bottomPG1T z_centerPG1T z_topPG3T z_bottomPG3T z_centerPG3T];
                    
                    t = (1:size(xnose, 1))'/FrameRate;
                    
                    if app.NoLightButton.Value
                        laserstart = [];
                        laserstop = [];
                    else
                        laser_timestamps = Bite_events(j).laser_timestamps;
                    end
                    
                    if ~isempty(laserstart)
                        if ~isempty(adjustment_start)
                            for k = 1:ninhibition
                                for m = 1:numel(adjustment_start)
                                    timestamp = adjustment_start(m);
                                    if timestamp >= laser_timestamps(2*k-1) && timestamp <= laser_timestamps(2*k)
                                        [~, id] = min(abs(t-timestamp));
                                        if id-trange >= 1
                                            ceyel = ceyelT(id, :);
%                                             ceyer = ceyerT(id, :);
                                            cnose = cnoseT(id, :);
%                                             if ~any(isnan([ceyel ceyer]))
                                            if ~any(isnan([ceyel cnose]))
%                                                 a = norm(ceyer);
                                                a = norm(cnose);
%                                                 b = norm(ceyel);
                                                b = norm(ceyel-cnose);
%                                                 c = norm(ceyel-ceyer);
                                                c = norm(ceyel);
%                                                 y = sqrt(a^2*c^2-dot(-ceyer, ceyel-ceyer)^2)/c; % new y coordinate of left & right eyes
                                                y = -sqrt(a^2*c^2-dot(cnose, ceyel)^2)/c; % new y coordinate of nose
%                                                 x1 = -dot(-ceyer, ceyel-ceyer)/c; % new x coordinate of right eye
                                                x1 = dot(cnose, ceyel)/c; % new x coordinate of nose
%                                                 x2 = dot(ceyel, ceyel-ceyer)/c; % new x coordinate of left eye
                                                x2 = c; % new x coordinate of left eye
                                                if 0
                                                    if x2 < 0.5 % x1 < -10 || x2 < 2.5 || x1 > -2.5 || x2 > 10 || x2-x1 < 5 || y < 11.5 || y > 13
                                                        continue;
                                                    end
                                                end
                                                if 0
                                                    if (i == 1 && j == 2 && m == 8) || (i == 1 && j == 2 && m == 21)
                                                        continue;
                                                    end
                                                end
                                                % standard othogonal bases of the new coordinate system
%                                                 u = (ceyel-ceyer)/c;
%                                                 v = (ceyer+(ceyel-ceyer)*(-x1)/c)/y;
%                                                 w = cross(ceyel, ceyer)/norm(cross(ceyel, ceyer));
                                                u = ceyel/c;
                                                v = (-cnose+ceyel/c*x1)/(-y);
                                                w = cross(cnose, ceyel)/norm(cross(cnose, ceyel));
                                                
                                                cpawl = cpawlT(id, :);
                                                cpawr = cpawrT(id, :);
                                                ctop = ctopT(id, :);
                                                cbottom = cbottomT(id, :);
                                                ccenter = ccenterT(id, :);
                                                % re-estimate pasta y coordinate based on z coordinate
                                                y_pasta = y_pastaT(id, :);
                                                z_pasta = z_pastaT(id, :);
                                                ctop = reestimate_ycoordinate_pasta(ctop, y_pasta, z_pasta);
                                                cbottom = reestimate_ycoordinate_pasta(cbottom, y_pasta, z_pasta);
                                                ccenter = reestimate_ycoordinate_pasta(ccenter, y_pasta, z_pasta);
                                                
                                                nadj_I(k) = nadj_I(k)+1;
                                                ceyel_all_I(:, :, nadj_I(k), k) = nan(2*trange+1, 3);
%                                                 ceyer_all_I(:, :, nadj_I(k), k) = nan(2*trange+1, 3);
                                                cnose_all_I(:, :, nadj_I(k), k) = nan(2*trange+1, 3);
                                                cpawl_all_I(:, :, nadj_I(k), k) = nan(2*trange+1, 3);
                                                cpawr_all_I(:, :, nadj_I(k), k) = nan(2*trange+1, 3);
                                                ctop_all_I(:, :, nadj_I(k), k) = nan(2*trange+1, 3);
                                                cbottom_all_I(:, :, nadj_I(k), k) = nan(2*trange+1, 3);
                                                ccenter_all_I(:, :, nadj_I(k), k) = nan(2*trange+1, 3);
                                                % compute the new coordinates
                                                ceyel_all_I(trange+1, :, nadj_I(k), k) = [u; v; w]*ceyel';
%                                                 ceyer_all_I(trange+1, :, nadj_I(k), k) = [u; v; w]*ceyer';
                                                cnose_all_I(trange+1, :, nadj_I(k), k) = [u; v; w]*cnose';
                                                cpawl_all_I(trange+1, :, nadj_I(k), k) = [u; v; w]*cpawl';
                                                cpawr_all_I(trange+1, :, nadj_I(k), k) = [u; v; w]*cpawr';
                                                ctop_all_I(trange+1, :, nadj_I(k), k) = [u; v; w]*ctop';
                                                cbottom_all_I(trange+1, :, nadj_I(k), k) = [u; v; w]*cbottom';
                                                ccenter_all_I(trange+1, :, nadj_I(k), k) = [u; v; w]*ccenter';
                                                
                                                adjinfo_I(nadj_I(k), k) = {['i = ' num2str(i) '; j = ' num2str(j) '; time = ' num2str(timestamp) '; m = ' num2str(m)]};
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    else
                        if ~isempty(adjustment_start)
                            for k = 1:ninhibition
                                for m = 1:numel(adjustment_start)
                                    timestamp = adjustment_start(m);
                                    if timestamp >= tinhibition(2*k-1) && timestamp <= tinhibition(2*k)
                                        [~, id] = min(abs(t-timestamp));
                                        if id-trange >= 1
                                            ceyel = ceyelT(id, :);
%                                             ceyer = ceyerT(id, :);
                                            cnose = cnoseT(id, :);
%                                             if ~any(isnan([ceyel ceyer]))
                                            if ~any(isnan([ceyel cnose]))
%                                                 a = norm(ceyer);
                                                a = norm(cnose);
%                                                 b = norm(ceyel);
                                                b = norm(ceyel-cnose);
%                                                 c = norm(ceyel-ceyer);
                                                c = norm(ceyel);
%                                                 y = sqrt(a^2*c^2-dot(-ceyer, ceyel-ceyer)^2)/c; % new y coordinate of left & right eyes
                                                y = -sqrt(a^2*c^2-dot(cnose, ceyel)^2)/c; % new y coordinate of nose
%                                                 x1 = -dot(-ceyer, ceyel-ceyer)/c; % new x coordinate of right eye
                                                x1 = dot(cnose, ceyel)/c; % new x coordinate of nose
%                                                 x2 = dot(ceyel, ceyel-ceyer)/c; % new x coordinate of left eye
                                                x2 = c; % new x coordinate of left eye
                                                if 0
                                                    if x2 < 0.5 % x1 < -10 || x2 < 2.5 || x1 > -2.5 || x2 > 10 || x2-x1 < 5 || y < 11.5 || y > 13
                                                        continue;
                                                    end
                                                end
                                                if 0
                                                    if (i == 1 && j == 2 && m == 8) || (i == 1 && j == 2 && m == 21)
                                                        continue;
                                                    end
                                                end
                                                % standard othogonal bases of the new coordinate system
%                                                 u = (ceyel-ceyer)/c;
%                                                 v = (ceyer+(ceyel-ceyer)*(-x1)/c)/y;
%                                                 w = cross(ceyel, ceyer)/norm(cross(ceyel, ceyer));
                                                u = ceyel/c;
                                                v = (-cnose+ceyel/c*x1)/(-y);
                                                w = cross(cnose, ceyel)/norm(cross(cnose, ceyel));
                                                
                                                cpawl = cpawlT(id, :);
                                                cpawr = cpawrT(id, :);
                                                ctop = ctopT(id, :);
                                                cbottom = cbottomT(id, :);
                                                ccenter = ccenterT(id, :);
                                                % re-estimate pasta y coordinate based on z coordinate
                                                y_pasta = y_pastaT(id, :);
                                                z_pasta = z_pastaT(id, :);
                                                ctop = reestimate_ycoordinate_pasta(ctop, y_pasta, z_pasta);
                                                cbottom = reestimate_ycoordinate_pasta(cbottom, y_pasta, z_pasta);
                                                ccenter = reestimate_ycoordinate_pasta(ccenter, y_pasta, z_pasta);
                                                
                                                nadj_NI(k) = nadj_NI(k)+1;
                                                ceyel_all_NI(:, :, nadj_NI(k), k) = nan(2*trange+1, 3);
%                                                 ceyer_all_NI(:, :, nadj_NI(k), k) = nan(2*trange+1, 3);
                                                cnose_all_NI(:, :, nadj_NI(k), k) = nan(2*trange+1, 3);
                                                cpawl_all_NI(:, :, nadj_NI(k), k) = nan(2*trange+1, 3);
                                                cpawr_all_NI(:, :, nadj_NI(k), k) = nan(2*trange+1, 3);
                                                ctop_all_NI(:, :, nadj_NI(k), k) = nan(2*trange+1, 3);
                                                cbottom_all_NI(:, :, nadj_NI(k), k) = nan(2*trange+1, 3);
                                                ccenter_all_NI(:, :, nadj_NI(k), k) = nan(2*trange+1, 3);
                                                % compute the new coordinates
                                                ceyel_all_NI(trange+1, :, nadj_NI(k), k) = [u; v; w]*ceyel';
%                                                 ceyer_all_NI(trange+1, :, nadj_NI(k), k) = [u; v; w]*ceyer';
                                                cnose_all_NI(trange+1, :, nadj_NI(k), k) = [u; v; w]*cnose';
                                                cpawl_all_NI(trange+1, :, nadj_NI(k), k) = [u; v; w]*cpawl';
                                                cpawr_all_NI(trange+1, :, nadj_NI(k), k) = [u; v; w]*cpawr';
                                                ctop_all_NI(trange+1, :, nadj_NI(k), k) = [u; v; w]*ctop';
                                                cbottom_all_NI(trange+1, :, nadj_NI(k), k) = [u; v; w]*cbottom';
                                                ccenter_all_NI(trange+1, :, nadj_NI(k), k) = [u; v; w]*ccenter';
                                                
                                                adjinfo_NI(nadj_NI(k), k) = {['i = ' num2str(i) '; j = ' num2str(j) '; time = ' num2str(timestamp) '; m = ' num2str(m)]};
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                else
                    [~, ~, ~, ~, ~, xanklel, yanklel, zanklel, ~, ~, ~, ~] = trajectory_postprocessing(17, Exp_Paths{i}, j, label_table, nmedian, FrameRate);
                    [~, ~, ~, ~, ~, xankler, yankler, zankler, ~, ~, ~, ~] = trajectory_postprocessing(18, Exp_Paths{i}, j, label_table, nmedian, FrameRate);
                    cnose = [];
                    cpawl = [];
                    cpawr = [];
                    ctop = [];
                    cbottom = [];
                    ccenter = [];
                    % new coordinate origin is average of left and right ankle
                    corigin = [(xanklel(:, 1)+xankler(:, 1))/2 (yanklel(:, 1)+yankler(:, 1))/2 (zanklel(:, 1)+zankler(:, 1))/2];
                    % v means vector
                    vankleRL = [xanklel(:, 1)-xankler(:, 1) yanklel(:, 1)-yankler(:, 1)];
                    theta = atan(vankleRL(:, 2)./vankleRL(:, 1));
                    theta(vankleRL(:, 1) < 0) = pi+theta(vankleRL(:, 1) < 0);
                    % coordinate transformation, translation first, then rotation
                    cnoseT = [xnose(:, 1)-corigin(:, 1) ynose(:, 1)-corigin(:, 2) znose(:, 1)-corigin(:, 3)];
                    cnose(:, 1) = cnoseT(:, 1).*cos(theta)+cnoseT(:, 2).*sin(theta);
                    cnose(:, 2) = cnoseT(:, 2).*cos(theta)-cnoseT(:, 1).*sin(theta);
                    cnose(:, 3) = cnoseT(:, 3);
                    cpawlT = [xpawl(:, 1)-corigin(:, 1) ypawl(:, 1)-corigin(:, 2) zpawl(:, 1)-corigin(:, 3)];
                    cpawl(:, 1) = cpawlT(:, 1).*cos(theta)+cpawlT(:, 2).*sin(theta);
                    cpawl(:, 2) = cpawlT(:, 2).*cos(theta)-cpawlT(:, 1).*sin(theta);
                    cpawl(:, 3) = cpawlT(:, 3);
                    cpawrT = [xpawr(:, 1)-corigin(:, 1) ypawr(:, 1)-corigin(:, 2) zpawr(:, 1)-corigin(:, 3)];
                    cpawr(:, 1) = cpawrT(:, 1).*cos(theta)+cpawrT(:, 2).*sin(theta);
                    cpawr(:, 2) = cpawrT(:, 2).*cos(theta)-cpawrT(:, 1).*sin(theta);
                    cpawr(:, 3) = cpawrT(:, 3);
                    
                    y_topPG1T = y_topPG1-corigin(:, 2);
                    y_topPG3T = y_topPG3-corigin(:, 2);
                    z_topPG1T = z_topPG1-corigin(:, 3);
                    z_topPG3T = z_topPG3-corigin(:, 3);
                    y_bottomPG1T = y_bottomPG1-corigin(:, 2);
                    y_bottomPG3T = y_bottomPG3-corigin(:, 2);
                    z_bottomPG1T = z_bottomPG1-corigin(:, 3);
                    z_bottomPG3T = z_bottomPG3-corigin(:, 3);
                    y_centerPG1T = y_centerPG1-corigin(:, 2);
                    y_centerPG3T = y_centerPG3-corigin(:, 2);
                    z_centerPG1T = z_centerPG1-corigin(:, 3);
                    z_centerPG3T = z_centerPG3-corigin(:, 3);
                    y_pastaT = [y_topPG1T y_bottomPG1T y_centerPG1T y_topPG3T y_bottomPG3T y_centerPG3T];
                    z_pastaT = [z_topPG1T z_bottomPG1T z_centerPG1T z_topPG3T z_bottomPG3T z_centerPG3T];
                    
                    ctopT = [xtop(:, 1)-corigin(:, 1) ytop(:, 1)-corigin(:, 2) ztop(:, 1)-corigin(:, 3)];
                    cbottomT = [xbottom(:, 1)-corigin(:, 1) ybottom(:, 1)-corigin(:, 2) zbottom(:, 1)-corigin(:, 3)];
                    ccenterT = [xcenter(:, 1)-corigin(:, 1) ycenter(:, 1)-corigin(:, 2) zcenter(:, 1)-corigin(:, 3)];
                    
                    t = (1:size(xnose, 1))'/FrameRate;
                    
                    if app.NoLightButton.Value
                        laserstart = [];
                        laserstop = [];
                    else
                        laser_timestamps = Bite_events(j).laser_timestamps;
                    end
                    
                    if ~isempty(laserstart)
                        if ~isempty(adjustment_start)
                            for k = 1:ninhibition
                                for m = 1:numel(adjustment_start)
                                    timestamp = adjustment_start(m);
                                    if timestamp >= laser_timestamps(2*k-1) && timestamp <= laser_timestamps(2*k)
                                        [~, id] = min(abs(t-timestamp));
                                        if id-trange >= 1
                                            ctopT_bite = ctopT(id, :);
                                            cbottomT_bite = cbottomT(id, :);
                                            ccenterT_bite = ccenterT(id, :);
                                            % re-estimate pasta y coordinate based on z coordinate
                                            y_pasta = y_pastaT(id, :);
                                            z_pasta = z_pastaT(id, :);
                                            ctopT_bite = reestimate_ycoordinate_pasta(ctopT_bite, y_pasta, z_pasta);
                                            cbottomT_bite = reestimate_ycoordinate_pasta(cbottomT_bite, y_pasta, z_pasta);
                                            ccenterT_bite = reestimate_ycoordinate_pasta(ccenterT_bite, y_pasta, z_pasta);

                                            ctop(1) = ctopT_bite(1)*cos(theta(id))+ctopT_bite(2)*sin(theta(id));
                                            ctop(2) = ctopT_bite(2)*cos(theta(id))-ctopT_bite(1)*sin(theta(id));
                                            ctop(3) = ctopT_bite(3);
                                            cbottom(1) = cbottomT_bite(1)*cos(theta(id))+cbottomT_bite(2)*sin(theta(id));
                                            cbottom(2) = cbottomT_bite(2)*cos(theta(id))-cbottomT_bite(1)*sin(theta(id));
                                            cbottom(3) = cbottomT_bite(3);
                                            ccenter(1) = ccenterT_bite(1)*cos(theta(id))+ccenterT_bite(2)*sin(theta(id));
                                            ccenter(2) = ccenterT_bite(2)*cos(theta(id))-ccenterT_bite(1)*sin(theta(id));
                                            ccenter(3) = ccenterT_bite(3);
                                            
                                            nadj_I(k) = nadj_I(k)+1;
                                            cnose_all_I(:, :, nadj_I(k), k) = nan(2*trange+1, 3);
                                            cpawl_all_I(:, :, nadj_I(k), k) = nan(2*trange+1, 3);
                                            cpawr_all_I(:, :, nadj_I(k), k) = nan(2*trange+1, 3);
                                            ctop_all_I(:, :, nadj_I(k), k) = nan(2*trange+1, 3);
                                            cbottom_all_I(:, :, nadj_I(k), k) = nan(2*trange+1, 3);
                                            ccenter_all_I(:, :, nadj_I(k), k) = nan(2*trange+1, 3);
                                            
                                            cnose_all_I(trange+1, :, nadj_I(k), k) = cnose(id, :);
                                            cpawl_all_I(trange+1, :, nadj_I(k), k) = cpawl(id, :);
                                            cpawr_all_I(trange+1, :, nadj_I(k), k) = cpawr(id, :);
                                            ctop_all_I(trange+1, :, nadj_I(k), k) = ctop;
                                            cbottom_all_I(trange+1, :, nadj_I(k), k) = cbottom;
                                            ccenter_all_I(trange+1, :, nadj_I(k), k) = ccenter;
                                        end
                                    end
                                end
                            end
                        end
                    else
                        if ~isempty(adjustment_start)
                            for k = 1:ninhibition
                                for m = 1:numel(adjustment_start)
                                    timestamp = adjustment_start(m);
                                    if timestamp >= tinhibition(2*k-1) && timestamp <= tinhibition(2*k)
                                        [~, id] = min(abs(t-timestamp));
                                        if id-trange >= 1
                                            ctopT_bite = ctopT(id, :);
                                            cbottomT_bite = cbottomT(id, :);
                                            ccenterT_bite = ccenterT(id, :);
                                            % re-estimate pasta y coordinate based on z coordinate
                                            y_pasta = y_pastaT(id, :);
                                            z_pasta = z_pastaT(id, :);
                                            ctopT_bite = reestimate_ycoordinate_pasta(ctopT_bite, y_pasta, z_pasta);
                                            cbottomT_bite = reestimate_ycoordinate_pasta(cbottomT_bite, y_pasta, z_pasta);
                                            ccenterT_bite = reestimate_ycoordinate_pasta(ccenterT_bite, y_pasta, z_pasta);
                                            
                                            ctop(1) = ctopT_bite(1)*cos(theta(id))+ctopT_bite(2)*sin(theta(id));
                                            ctop(2) = ctopT_bite(2)*cos(theta(id))-ctopT_bite(1)*sin(theta(id));
                                            ctop(3) = ctopT_bite(3);
                                            cbottom(1) = cbottomT_bite(1)*cos(theta(id))+cbottomT_bite(2)*sin(theta(id));
                                            cbottom(2) = cbottomT_bite(2)*cos(theta(id))-cbottomT_bite(1)*sin(theta(id));
                                            cbottom(3) = cbottomT_bite(3);
                                            ccenter(1) = ccenterT_bite(1)*cos(theta(id))+ccenterT_bite(2)*sin(theta(id));
                                            ccenter(2) = ccenterT_bite(2)*cos(theta(id))-ccenterT_bite(1)*sin(theta(id));
                                            ccenter(3) = ccenterT_bite(3);
                                            
                                            nadj_NI(k) = nadj_NI(k)+1;
                                            cnose_all_NI(:, :, nadj_NI(k), k) = nan(2*trange+1, 3);
                                            cpawl_all_NI(:, :, nadj_NI(k), k) = nan(2*trange+1, 3);
                                            cpawr_all_NI(:, :, nadj_NI(k), k) = nan(2*trange+1, 3);
                                            ctop_all_NI(:, :, nadj_NI(k), k) = nan(2*trange+1, 3);
                                            cbottom_all_NI(:, :, nadj_NI(k), k) = nan(2*trange+1, 3);
                                            ccenter_all_NI(:, :, nadj_NI(k), k) = nan(2*trange+1, 3);
                                            
                                            cnose_all_NI(trange+1, :, nadj_NI(k), k) = cnose(id, :);
                                            cpawl_all_NI(trange+1, :, nadj_NI(k), k) = cpawl(id, :);
                                            cpawr_all_NI(trange+1, :, nadj_NI(k), k) = cpawr(id, :);
                                            ctop_all_NI(trange+1, :, nadj_NI(k), k) = ctop;
                                            cbottom_all_NI(trange+1, :, nadj_NI(k), k) = cbottom;
                                            ccenter_all_NI(trange+1, :, nadj_NI(k), k) = ccenter;
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end
if analyze_method == 2
    t = (-trange:1:trange)'/FrameRate;
    if ~app.NoLightButton.Value
        for i = 1:ninhibition
            figure;
            plot_tj_acrosstrial(t, znoseadj_NI(:, :, i), [0 0 0], [0 0 0], 'Z (mm)');
            plot_tj_acrosstrial(t, znoseadj_I(:, :, i), [0.5 0.5 0.5], [0.5 0.5 0.5], 'Z (mm)');
            plot_tj_acrosstrial(t, zpawladj_NI(:, :, i), [0 0 0.5], [0 0 0.5], 'Z (mm)');
            plot_tj_acrosstrial(t, zpawladj_I(:, :, i), [0 0 1], [0 0 1], 'Z (mm)');
            plot_tj_acrosstrial(t, zpawradj_NI(:, :, i), [0 0.5 0], [0 0.5 0], 'Z (mm)');
            plot_tj_acrosstrial(t, zpawradj_I(:, :, i), [0 1 0], [0 1 0], 'Z (mm)');
            yl = ylim;
            plot([0 0], yl, '--k');
            
            figure;
            plot_tj_acrosstrial(t, znoseadj_NI(:, :, i)-zpawladj_NI(:, :, i), [0 0 0.5], [0 0 0.5], '\Delta Z (mm)');
            plot_tj_acrosstrial(t, znoseadj_I(:, :, i)-zpawladj_I(:, :, i), [0 0 1], [0 0 1], '\Delta Z (mm)');
            plot_tj_acrosstrial(t, znoseadj_NI(:, :, i)-zpawradj_NI(:, :, i), [0 0.5 0], [0 0.5 0], '\Delta Z (mm)');
            plot_tj_acrosstrial(t, znoseadj_I(:, :, i)-zpawradj_I(:, :, i), [0 1 0], [0 1 0], '\Delta Z (mm)');
            plot_tj_acrosstrial(t, zpawladj_NI(:, :, i)-zpawradj_NI(:, :, i), [0 0 0], [0 0 0], '\Delta Z (mm)');
            plot_tj_acrosstrial(t, zpawladj_I(:, :, i)-zpawradj_I(:, :, i), [0.5 0.5 0.5], [0.5 0.5 0.5], '\Delta Z (mm)');
            yl = ylim;
            plot([0 0], yl, '--k');
            
            figure;
            subplot(1, 3, 1);
            plot_tj_acrosstrial(t, zspeednoseadj_NI(:, :, i), [0 0 0], [0 0 0], 'Z speed (mm/s)');
            plot_tj_acrosstrial(t, zspeednoseadj_I(:, :, i), [0.5 0.5 0.5], [0.5 0.5 0.5], 'Z speed (mm/s)');
            yl = ylim;
            plot([0 0], yl, '--k');
            subplot(1, 3, 2);
            plot_tj_acrosstrial(t, zspeedpawladj_NI(:, :, i), [0 0 0.5], [0 0 0.5], 'Z speed (mm/s)');
            plot_tj_acrosstrial(t, zspeedpawladj_I(:, :, i), [0 0 1], [0 0 1], 'Z speed (mm/s)');
            yl = ylim;
            plot([0 0], yl, '--k');
            subplot(1, 3, 3);
            plot_tj_acrosstrial(t, zspeedpawradj_NI(:, :, i), [0 0.5 0], [0 0.5 0], 'Z speed (mm/s)');
            plot_tj_acrosstrial(t, zspeedpawradj_I(:, :, i), [0 1 0], [0 1 0], 'Z speed (mm/s)');
            yl = ylim;
            plot([0 0], yl, '--k');
            
            figure;
            subplot(1, 3, 1);
            plot_tj_acrosstrial(t, xyzspeednoseadj_NI(:, :, i), [0 0 0], [0 0 0], 'XYZ speed (mm/s)');
            plot_tj_acrosstrial(t, xyzspeednoseadj_I(:, :, i), [0.5 0.5 0.5], [0.5 0.5 0.5], 'XYZ speed (mm/s)');
            yl = ylim;
            plot([0 0], yl, '--k');
            subplot(1, 3, 2);
            plot_tj_acrosstrial(t, xyzspeedpawladj_NI(:, :, i), [0 0 0.5], [0 0 0.5], 'XYZ speed (mm/s)');
            plot_tj_acrosstrial(t, xyzspeedpawladj_I(:, :, i), [0 0 1], [0 0 1], 'XYZ speed (mm/s)');
            yl = ylim;
            plot([0 0], yl, '--k');
            subplot(1, 3, 3);
            plot_tj_acrosstrial(t, xyzspeedpawradj_NI(:, :, i), [0 0.5 0], [0 0.5 0], 'XYZ speed (mm/s)');
            plot_tj_acrosstrial(t, xyzspeedpawradj_I(:, :, i), [0 1 0], [0 1 0], 'XYZ speed (mm/s)');
            yl = ylim;
            plot([0 0], yl, '--k');
        end
    else
        for i = 1:ninhibition
            figure;
            plot_tj_acrosstrial(t, znoseadj_NI(:, :, i), [0 0 0], [0 0 0], 'Z (mm)');
            plot_tj_acrosstrial(t, zpawladj_NI(:, :, i), [0 0 0.5], [0 0 0.5], 'Z (mm)');
            plot_tj_acrosstrial(t, zpawradj_NI(:, :, i), [0 0.5 0], [0 0.5 0], 'Z (mm)');
            yl = ylim;
            plot([0 0], yl, '--k');
            
            figure;
            plot_tj_acrosstrial(t, znoseadj_NI(:, :, i)-zpawladj_NI(:, :, i), [0 0 0.5], [0 0 0.5], '\Delta Z (mm)');
            plot_tj_acrosstrial(t, znoseadj_NI(:, :, i)-zpawradj_NI(:, :, i), [0 0.5 0], [0 0.5 0], '\Delta Z (mm)');
            plot_tj_acrosstrial(t, zpawladj_NI(:, :, i)-zpawradj_NI(:, :, i), [0 0 0], [0 0 0], '\Delta Z (mm)');
            yl = ylim;
            plot([0 0], yl, '--k');
            
            figure;
            subplot(1, 3, 1);
            plot_tj_acrosstrial(t, zspeednoseadj_NI(:, :, i), [0 0 0], [0 0 0], 'Z speed (mm/s)');
            yl = ylim;
            plot([0 0], yl, '--k');
            subplot(1, 3, 2);
            plot_tj_acrosstrial(t, zspeedpawladj_NI(:, :, i), [0 0 0.5], [0 0 0.5], 'Z speed (mm/s)');
            yl = ylim;
            plot([0 0], yl, '--k');
            subplot(1, 3, 3);
            plot_tj_acrosstrial(t, zspeedpawradj_NI(:, :, i), [0 0.5 0], [0 0.5 0], 'Z speed (mm/s)');
            yl = ylim;
            plot([0 0], yl, '--k');
            
            figure;
            subplot(1, 3, 1);
            plot_tj_acrosstrial(t, xyzspeednoseadj_NI(:, :, i), [0 0 0], [0 0 0], 'XYZ speed (mm/s)');
            yl = ylim;
            plot([0 0], yl, '--k');
            subplot(1, 3, 2);
            plot_tj_acrosstrial(t, xyzspeedpawladj_NI(:, :, i), [0 0 0.5], [0 0 0.5], 'XYZ speed (mm/s)');
            yl = ylim;
            plot([0 0], yl, '--k');
            subplot(1, 3, 3);
            plot_tj_acrosstrial(t, xyzspeedpawradj_NI(:, :, i), [0 0.5 0], [0 0.5 0], 'XYZ speed (mm/s)');
            yl = ylim;
            plot([0 0], yl, '--k');
        end
    end
elseif analyze_method == 1 || analyze_method == 3
    if analyze_method == 1
        if ~app.NoLightButton.Value
            for i = 1:ninhibition
                % separate the pasta from the two different sides
                if nadj_NI(i) ~= 0
                    [IDleft_front_NI, IDright_front_NI, IDleft_top_NI, IDright_top_NI] = separate_pasta_sides(ctop_all_NI(trange+1, :, :, i), cbottom_all_NI(trange+1, :, :, i));
                else
                    IDleft_front_NI = [];
                    IDright_front_NI = [];
                    IDleft_top_NI = [];
                    IDright_top_NI = [];
                end
                if nadj_I(i) ~= 0
                    [IDleft_front_I, IDright_front_I, IDleft_top_I, IDright_top_I] = separate_pasta_sides(ctop_all_I(trange+1, :, :, i), cbottom_all_I(trange+1, :, :, i));
                else
                    IDleft_front_I = [];
                    IDright_front_I = [];
                    IDleft_top_I = [];
                    IDright_top_I = [];
                end
                
%                 draw_individual_locations_at_bite(ceyel_all_NI, ceyer_all_NI, cpawl_all_NI, cpawr_all_NI, ctop_all_NI, cbottom_all_NI,...
%                     ceyel_all_I, ceyer_all_I, cpawl_all_I, cpawr_all_I, ctop_all_I, cbottom_all_I, [IDleft_front_NI IDright_front_NI], [IDleft_top_NI IDright_top_NI],...
%                     [IDleft_front_I IDright_front_I], [IDleft_top_I IDright_top_I], trange, i, adjinfo_NI, adjinfo_I);
                draw_individual_locations_at_bite_new(ceyel_all_NI, cnose_all_NI, cpawl_all_NI, cpawr_all_NI, ctop_all_NI, cbottom_all_NI,...
                    ceyel_all_I, cnose_all_I, cpawl_all_I, cpawr_all_I, ctop_all_I, cbottom_all_I, [IDleft_front_NI IDright_front_NI], [IDleft_top_NI IDright_top_NI],...
                    [IDleft_front_I IDright_front_I], [IDleft_top_I IDright_top_I], trange, i, adjinfo_NI, adjinfo_I);
                set(gcf, 'Name', 'All');
%                 draw_individual_locations_at_bite(ceyel_all_NI, ceyer_all_NI, cpawl_all_NI, cpawr_all_NI, ctop_all_NI, cbottom_all_NI,...
%                     ceyel_all_I, ceyer_all_I, cpawl_all_I, cpawr_all_I, ctop_all_I, cbottom_all_I, IDleft_front_NI, IDleft_top_NI, IDleft_front_I, IDleft_top_I, trange, i, adjinfo_NI, adjinfo_I);
                draw_individual_locations_at_bite_new(ceyel_all_NI, cnose_all_NI, cpawl_all_NI, cpawr_all_NI, ctop_all_NI, cbottom_all_NI,...
                    ceyel_all_I, cnose_all_I, cpawl_all_I, cpawr_all_I, ctop_all_I, cbottom_all_I, IDleft_front_NI, IDleft_top_NI, IDleft_front_I, IDleft_top_I, trange, i, adjinfo_NI, adjinfo_I);
                set(gcf, 'Name', 'Left');
%                 draw_individual_locations_at_bite(ceyel_all_NI, ceyer_all_NI, cpawl_all_NI, cpawr_all_NI, ctop_all_NI, cbottom_all_NI,...
%                     ceyel_all_I, ceyer_all_I, cpawl_all_I, cpawr_all_I, ctop_all_I, cbottom_all_I, IDright_front_NI, IDright_top_NI, IDright_front_I, IDright_top_I, trange, i, adjinfo_NI, adjinfo_I);
                draw_individual_locations_at_bite_new(ceyel_all_NI, cnose_all_NI, cpawl_all_NI, cpawr_all_NI, ctop_all_NI, cbottom_all_NI,...
                    ceyel_all_I, cnose_all_I, cpawl_all_I, cpawr_all_I, ctop_all_I, cbottom_all_I, IDright_front_NI, IDright_top_NI, IDright_front_I, IDright_top_I, trange, i, adjinfo_NI, adjinfo_I);
                set(gcf, 'Name', 'Right');
                
                % mean and errorbar
%                 draw_average_locations_at_bite(ceyel_all_NI, ceyer_all_NI, cpawl_all_NI, cpawr_all_NI, ctop_all_NI, cbottom_all_NI,...
%                     ceyel_all_I, ceyer_all_I, cpawl_all_I, cpawr_all_I, ctop_all_I, cbottom_all_I, [IDleft_front_NI IDright_front_NI], [IDleft_top_NI IDright_top_NI],...
%                     [IDleft_front_I IDright_front_I], [IDleft_top_I IDright_top_I], trange, i);
                draw_average_locations_at_bite_new(ceyel_all_NI, cnose_all_NI, cpawl_all_NI, cpawr_all_NI, ctop_all_NI, cbottom_all_NI,...
                    ceyel_all_I, cnose_all_I, cpawl_all_I, cpawr_all_I, ctop_all_I, cbottom_all_I, [IDleft_front_NI IDright_front_NI], [IDleft_top_NI IDright_top_NI],...
                    [IDleft_front_I IDright_front_I], [IDleft_top_I IDright_top_I], trange, i);
                set(gcf, 'Name', 'All');

%                 draw_average_locations_at_bite(ceyel_all_NI, ceyer_all_NI, cpawl_all_NI, cpawr_all_NI, ctop_all_NI, cbottom_all_NI,...
%                     ceyel_all_I, ceyer_all_I, cpawl_all_I, cpawr_all_I, ctop_all_I, cbottom_all_I, IDleft_front_NI, IDleft_top_NI, IDleft_front_I, IDleft_top_I, trange, i);
                draw_average_locations_at_bite_new(ceyel_all_NI, cnose_all_NI, cpawl_all_NI, cpawr_all_NI, ctop_all_NI, cbottom_all_NI,...
                    ceyel_all_I, cnose_all_I, cpawl_all_I, cpawr_all_I, ctop_all_I, cbottom_all_I, IDleft_front_NI, IDleft_top_NI, IDleft_front_I, IDleft_top_I, trange, i);
                set(gcf, 'Name', 'Left');

%                 draw_average_locations_at_bite(ceyel_all_NI, ceyer_all_NI, cpawl_all_NI, cpawr_all_NI, ctop_all_NI, cbottom_all_NI,...
%                     ceyel_all_I, ceyer_all_I, cpawl_all_I, cpawr_all_I, ctop_all_I, cbottom_all_I, IDright_front_NI, IDright_top_NI, IDright_front_I, IDright_top_I, trange, i);
                draw_average_locations_at_bite_new(ceyel_all_NI, cnose_all_NI, cpawl_all_NI, cpawr_all_NI, ctop_all_NI, cbottom_all_NI,...
                    ceyel_all_I, cnose_all_I, cpawl_all_I, cpawr_all_I, ctop_all_I, cbottom_all_I, IDright_front_NI, IDright_top_NI, IDright_front_I, IDright_top_I, trange, i);
                set(gcf, 'Name', 'Right');
                
                [orientation_NI, orientation_I] = draw_individual_pasta_at_bite(ctop_all_NI, cbottom_all_NI, ctop_all_I, cbottom_all_I,...
                    [IDleft_front_NI IDright_front_NI], [IDleft_front_I IDright_front_I], trange, i, adjinfo_NI, adjinfo_I);
                
                result.ceyel_all_NI = ceyel_all_NI;
%                 result.ceyer_all_NI = ceyer_all_NI;
                result.cnose_all_NI = cnose_all_NI;
                result.cpawl_all_NI = cpawl_all_NI;
                result.cpawr_all_NI = cpawr_all_NI;
                result.ctop_all_NI = ctop_all_NI;
                result.cbottom_all_NI = cbottom_all_NI;
                result.IDleft_front_NI = IDleft_front_NI;
                result.IDleft_top_NI = IDleft_top_NI;
                result.IDright_front_NI = IDright_front_NI;
                result.IDright_top_NI = IDright_top_NI;
                result.orientation_NI = orientation_NI;
                
                result.ceyel_all_I = ceyel_all_I;
%                 result.ceyer_all_I = ceyer_all_I;
                result.cnose_all_I = cnose_all_I;
                result.cpawl_all_I = cpawl_all_I;
                result.cpawr_all_I = cpawr_all_I;
                result.ctop_all_I = ctop_all_I;
                result.cbottom_all_I = cbottom_all_I;
                result.IDleft_front_I = IDleft_front_I;
                result.IDleft_top_I = IDleft_top_I;
                result.IDright_front_I = IDright_front_I;
                result.IDright_top_I = IDright_top_I;
                result.orientation_I = orientation_I;
                
                curpwd = pwd;
                try
                    cd(path);
                end
                [file, path] = uiputfile('AdjustmentState.mat', 'Save result');
                if file ~= 0
                    save([path file], 'result');
                end
                cd(curpwd);
            end
        else
            for i = 1:ninhibition
                % separate the pasta from the two different sides
                if nadj_NI(i) ~= 0
                    [IDleft_front_NI, IDright_front_NI, IDleft_top_NI, IDright_top_NI] = separate_pasta_sides(ctop_all_NI(trange+1, :, :, i), cbottom_all_NI(trange+1, :, :, i));
                else
                    IDleft_front_NI = [];
                    IDright_front_NI = [];
                    IDleft_top_NI = [];
                    IDright_top_NI = [];
                end
                
%                 draw_individual_locations_at_bite(ceyel_all_NI, ceyer_all_NI, cpawl_all_NI, cpawr_all_NI, ctop_all_NI, cbottom_all_NI,...
%                     [], [], [], [], [], [], [IDleft_front_NI IDright_front_NI], [IDleft_top_NI IDright_top_NI], [], [], trange, i, adjinfo_NI, []);
                draw_individual_locations_at_bite_new(ceyel_all_NI, cnose_all_NI, cpawl_all_NI, cpawr_all_NI, ctop_all_NI, cbottom_all_NI,...
                    [], [], [], [], [], [], [IDleft_front_NI IDright_front_NI], [IDleft_top_NI IDright_top_NI], [], [], trange, i, adjinfo_NI, []);
                set(gcf, 'Name', 'All');
%                 draw_individual_locations_at_bite(ceyel_all_NI, ceyer_all_NI, cpawl_all_NI, cpawr_all_NI, ctop_all_NI, cbottom_all_NI,...
%                     [], [], [], [], [], [], IDleft_front_NI, IDleft_top_NI, [], [], trange, i, adjinfo_NI, []);
                draw_individual_locations_at_bite_new(ceyel_all_NI, cnose_all_NI, cpawl_all_NI, cpawr_all_NI, ctop_all_NI, cbottom_all_NI,...
                    [], [], [], [], [], [], IDleft_front_NI, IDleft_top_NI, [], [], trange, i, adjinfo_NI, []);
                set(gcf, 'Name', 'Left');
%                 draw_individual_locations_at_bite(ceyel_all_NI, ceyer_all_NI, cpawl_all_NI, cpawr_all_NI, ctop_all_NI, cbottom_all_NI,...
%                     [], [], [], [], [], [], IDright_front_NI, IDright_top_NI, [], [], trange, i, adjinfo_NI, []);
                draw_individual_locations_at_bite_new(ceyel_all_NI, cnose_all_NI, cpawl_all_NI, cpawr_all_NI, ctop_all_NI, cbottom_all_NI,...
                    [], [], [], [], [], [], IDright_front_NI, IDright_top_NI, [], [], trange, i, adjinfo_NI, []);
                set(gcf, 'Name', 'Right');
                
                % mean and errorbar
%                 draw_average_locations_at_bite(ceyel_all_NI, ceyer_all_NI, cpawl_all_NI, cpawr_all_NI, ctop_all_NI, cbottom_all_NI,...
%                     [], [], [], [], [], [], [IDleft_front_NI IDright_front_NI], [IDleft_top_NI IDright_top_NI], [], [], trange, i);
                draw_average_locations_at_bite_new(ceyel_all_NI, cnose_all_NI, cpawl_all_NI, cpawr_all_NI, ctop_all_NI, cbottom_all_NI,...
                    [], [], [], [], [], [], [IDleft_front_NI IDright_front_NI], [IDleft_top_NI IDright_top_NI], [], [], trange, i);
                set(gcf, 'Name', 'All');

%                 draw_average_locations_at_bite(ceyel_all_NI, ceyer_all_NI, cpawl_all_NI, cpawr_all_NI, ctop_all_NI, cbottom_all_NI,...
%                     [], [], [], [], [], [], IDleft_front_NI, IDleft_top_NI, [], [], trange, i);
                draw_average_locations_at_bite_new(ceyel_all_NI, cnose_all_NI, cpawl_all_NI, cpawr_all_NI, ctop_all_NI, cbottom_all_NI,...
                    [], [], [], [], [], [], IDleft_front_NI, IDleft_top_NI, [], [], trange, i);
                set(gcf, 'Name', 'Left');

%                 draw_average_locations_at_bite(ceyel_all_NI, ceyer_all_NI, cpawl_all_NI, cpawr_all_NI, ctop_all_NI, cbottom_all_NI,...
%                     [], [], [], [], [], [], IDright_front_NI, IDright_top_NI, [], [], trange, i);
                draw_average_locations_at_bite_new(ceyel_all_NI, cnose_all_NI, cpawl_all_NI, cpawr_all_NI, ctop_all_NI, cbottom_all_NI,...
                    [], [], [], [], [], [], IDright_front_NI, IDright_top_NI, [], [], trange, i);
                set(gcf, 'Name', 'Right');

                [orientation_NI, ~] = draw_individual_pasta_at_bite(ctop_all_NI, cbottom_all_NI, [], [], [IDleft_front_NI IDright_front_NI], [], trange, i, adjinfo_NI, []);
                
                result.ceyel_all_NI = ceyel_all_NI;
%                 result.ceyer_all_NI = ceyer_all_NI;
                result.cnose_all_NI = cnose_all_NI;
                result.cpawl_all_NI = cpawl_all_NI;
                result.cpawr_all_NI = cpawr_all_NI;
                result.ctop_all_NI = ctop_all_NI;
                result.cbottom_all_NI = cbottom_all_NI;
                result.IDleft_front_NI = IDleft_front_NI;
                result.IDleft_top_NI = IDleft_top_NI;
                result.IDright_front_NI = IDright_front_NI;
                result.IDright_top_NI = IDright_top_NI;
                result.orientation_NI = orientation_NI;
                
                curpwd = pwd;
                try
                    cd(path);
                end
                [file, path] = uiputfile('AdjustmentState.mat', 'Save result');
                if file ~= 0
                    save([path file], 'result');
                end
                cd(curpwd);
            end
        end
    else
        if ~app.NoLightButton.Value
            for i = 1:ninhibition
                figure;
                subplot(1, 2, 1);
                hold on;
                for j = 1:nadj_NI(i)
                    if ~any(isnan([cnose_all_NI(trange+1, [1 3], j, i) cpawl_all_NI(trange+1, [1 3], j, i) cpawr_all_NI(trange+1, [1 3], j, i)]))
                        plot([cnose_all_NI(trange+1, 1, j, i) cpawl_all_NI(trange+1, 1, j, i)], [cnose_all_NI(trange+1, 3, j, i) cpawl_all_NI(trange+1, 3, j, i)], '-k');
                        plot([cnose_all_NI(trange+1, 1, j, i) cpawr_all_NI(trange+1, 1, j, i)], [cnose_all_NI(trange+1, 3, j, i) cpawr_all_NI(trange+1, 3, j, i)], '-k');
                        plot([cpawl_all_NI(trange+1, 1, j, i) cpawr_all_NI(trange+1, 1, j, i)], [cpawl_all_NI(trange+1, 3, j, i) cpawr_all_NI(trange+1, 3, j, i)], '-k');
                        plot(cnose_all_NI(trange+1, 1, j, i), cnose_all_NI(trange+1, 3, j, i), '.b', 'MarkerSize', 12);
                        plot(cpawl_all_NI(trange+1, 1, j, i), cpawl_all_NI(trange+1, 3, j, i), '.c', 'MarkerSize', 12);
                        plot(cpawr_all_NI(trange+1, 1, j, i), cpawr_all_NI(trange+1, 3, j, i), '.g', 'MarkerSize', 12);
                    end
                    if ~any(isnan([ctop_all_NI(trange+1, [1 3], j, i) cbottom_all_NI(trange+1, [1 3], j, i)]))
                        plot([ctop_all_NI(trange+1, 1, j, i) cbottom_all_NI(trange+1, 1, j, i)], [ctop_all_NI(trange+1, 3, j, i) cbottom_all_NI(trange+1, 3, j, i)], '-k');
                        plot(ctop_all_NI(trange+1, 1, j, i), ctop_all_NI(trange+1, 3, j, i), '.', 'Color', [1 1 0], 'MarkerSize', 12);
                        plot(cbottom_all_NI(trange+1, 1, j, i), cbottom_all_NI(trange+1, 3, j, i), '.', 'Color', [1 1 0], 'MarkerSize', 12);
                    end
                end
                for j = 1:nadj_I(i)
                    if ~any(isnan([cnose_all_I(trange+1, [1 3], j, i) cpawl_all_I(trange+1, [1 3], j, i) cpawr_all_I(trange+1, [1 3], j, i)]))
                        plot([cnose_all_I(trange+1, 1, j, i) cpawl_all_I(trange+1, 1, j, i)], [cnose_all_I(trange+1, 3, j, i) cpawl_all_I(trange+1, 3, j, i)], '-b');
                        plot([cnose_all_I(trange+1, 1, j, i) cpawr_all_I(trange+1, 1, j, i)], [cnose_all_I(trange+1, 3, j, i) cpawr_all_I(trange+1, 3, j, i)], '-b');
                        plot([cpawl_all_I(trange+1, 1, j, i) cpawr_all_I(trange+1, 1, j, i)], [cpawl_all_I(trange+1, 3, j, i) cpawr_all_I(trange+1, 3, j, i)], '-b');
                        plot(cnose_all_I(trange+1, 1, j, i), cnose_all_I(trange+1, 3, j, i), '.b', 'MarkerSize', 12);
                        plot(cpawl_all_I(trange+1, 1, j, i), cpawl_all_I(trange+1, 3, j, i), '.c', 'MarkerSize', 12);
                        plot(cpawr_all_I(trange+1, 1, j, i), cpawr_all_I(trange+1, 3, j, i), '.g', 'MarkerSize', 12);
                    end
                    if ~any(isnan([ctop_all_I(trange+1, [1 3], j, i) cbottom_all_I(trange+1, [1 3], j, i)]))
                        plot([ctop_all_I(trange+1, 1, j, i) cbottom_all_I(trange+1, 1, j, i)], [ctop_all_I(trange+1, 3, j, i) cbottom_all_I(trange+1, 3, j, i)], '-b');
                        plot(ctop_all_I(trange+1, 1, j, i), ctop_all_I(trange+1, 3, j, i), '.', 'Color', [1 1 0], 'MarkerSize', 12);
                        plot(cbottom_all_I(trange+1, 1, j, i), cbottom_all_I(trange+1, 3, j, i), '.', 'Color', [1 1 0], 'MarkerSize', 12);
                    end
                end
                axis equal;
                xlabel('X (mm)');
                ylabel('Z (mm)');
                set(gca, 'TickDir', 'out', 'FontSize', 12, 'ButtonDownFcn', @extract_figure);
                
                subplot(1, 2, 2);
                hold on;
                for j = 1:nadj_NI(i)
                    if ~any(isnan([cnose_all_NI(trange+1, [1 2], j, i) cpawl_all_NI(trange+1, [1 2], j, i) cpawr_all_NI(trange+1, [1 2], j, i)]))
                        plot([cnose_all_NI(trange+1, 1, j, i) cpawl_all_NI(trange+1, 1, j, i)], [cnose_all_NI(trange+1, 2, j, i) cpawl_all_NI(trange+1, 2, j, i)], '-k');
                        plot([cnose_all_NI(trange+1, 1, j, i) cpawr_all_NI(trange+1, 1, j, i)], [cnose_all_NI(trange+1, 2, j, i) cpawr_all_NI(trange+1, 2, j, i)], '-k');
                        plot([cpawl_all_NI(trange+1, 1, j, i) cpawr_all_NI(trange+1, 1, j, i)], [cpawl_all_NI(trange+1, 2, j, i) cpawr_all_NI(trange+1, 2, j, i)], '-k');
                        plot(cnose_all_NI(trange+1, 1, j, i), cnose_all_NI(trange+1, 2, j, i), '.b', 'MarkerSize', 12);
                        plot(cpawl_all_NI(trange+1, 1, j, i), cpawl_all_NI(trange+1, 2, j, i), '.c', 'MarkerSize', 12);
                        plot(cpawr_all_NI(trange+1, 1, j, i), cpawr_all_NI(trange+1, 2, j, i), '.g', 'MarkerSize', 12);
                    end
                    if ~any(isnan([ctop_all_NI(trange+1, [1 2], j, i) cbottom_all_NI(trange+1, [1 2], j, i)]))
                        plot([ctop_all_NI(trange+1, 1, j, i) cbottom_all_NI(trange+1, 1, j, i)], [ctop_all_NI(trange+1, 2, j, i) cbottom_all_NI(trange+1, 2, j, i)], '-k');
                        plot(ctop_all_NI(trange+1, 1, j, i), ctop_all_NI(trange+1, 2, j, i), '.', 'Color', [1 1 0], 'MarkerSize', 12);
                        plot(cbottom_all_NI(trange+1, 1, j, i), cbottom_all_NI(trange+1, 2, j, i), '.', 'Color', [1 1 0], 'MarkerSize', 12);
                    end
                end
                for j = 1:nadj_I(i)
                    if ~any(isnan([cnose_all_I(trange+1, [1 2], j, i) cpawl_all_I(trange+1, [1 2], j, i) cpawr_all_I(trange+1, [1 2], j, i)]))
                        plot([cnose_all_I(trange+1, 1, j, i) cpawl_all_I(trange+1, 1, j, i)], [cnose_all_I(trange+1, 2, j, i) cpawl_all_I(trange+1, 2, j, i)], '-b');
                        plot([cnose_all_I(trange+1, 1, j, i) cpawr_all_I(trange+1, 1, j, i)], [cnose_all_I(trange+1, 2, j, i) cpawr_all_I(trange+1, 2, j, i)], '-b');
                        plot([cpawl_all_I(trange+1, 1, j, i) cpawr_all_I(trange+1, 1, j, i)], [cpawl_all_I(trange+1, 2, j, i) cpawr_all_I(trange+1, 2, j, i)], '-b');
                        plot(cnose_all_I(trange+1, 1, j, i), cnose_all_I(trange+1, 2, j, i), '.b', 'MarkerSize', 12);
                        plot(cpawl_all_I(trange+1, 1, j, i), cpawl_all_I(trange+1, 2, j, i), '.c', 'MarkerSize', 12);
                        plot(cpawr_all_I(trange+1, 1, j, i), cpawr_all_I(trange+1, 2, j, i), '.g', 'MarkerSize', 12);
                    end
                    if ~any(isnan([ctop_all_I(trange+1, [1 2], j, i) cbottom_all_I(trange+1, [1 2], j, i)]))
                        plot([ctop_all_I(trange+1, 1, j, i) cbottom_all_I(trange+1, 1, j, i)], [ctop_all_I(trange+1, 2, j, i) cbottom_all_I(trange+1, 2, j, i)], '-b');
                        plot(ctop_all_I(trange+1, 1, j, i), ctop_all_I(trange+1, 2, j, i), '.', 'Color', [1 1 0], 'MarkerSize', 12);
                        plot(cbottom_all_I(trange+1, 1, j, i), cbottom_all_I(trange+1, 2, j, i), '.', 'Color', [1 1 0], 'MarkerSize', 12);
                    end
                end
                axis equal;
                xlabel('X (mm)');
                ylabel('Y (mm)');
                set(gca, 'TickDir', 'out', 'FontSize', 12, 'ButtonDownFcn', @extract_figure);
            end
        else
            for i = 1:ninhibition
                figure;
                subplot(1, 2, 1);
                hold on;
                for j = 1:nadj_NI(i)
                    if ~any(isnan([cnose_all_NI(trange+1, [1 3], j, i) cpawl_all_NI(trange+1, [1 3], j, i) cpawr_all_NI(trange+1, [1 3], j, i)]))
                        plot([cnose_all_NI(trange+1, 1, j, i) cpawl_all_NI(trange+1, 1, j, i)], [cnose_all_NI(trange+1, 3, j, i) cpawl_all_NI(trange+1, 3, j, i)], '-k');
                        plot([cnose_all_NI(trange+1, 1, j, i) cpawr_all_NI(trange+1, 1, j, i)], [cnose_all_NI(trange+1, 3, j, i) cpawr_all_NI(trange+1, 3, j, i)], '-k');
                        plot([cpawl_all_NI(trange+1, 1, j, i) cpawr_all_NI(trange+1, 1, j, i)], [cpawl_all_NI(trange+1, 3, j, i) cpawr_all_NI(trange+1, 3, j, i)], '-k');
                        plot(cnose_all_NI(trange+1, 1, j, i), cnose_all_NI(trange+1, 3, j, i), '.b', 'MarkerSize', 12);
                        plot(cpawl_all_NI(trange+1, 1, j, i), cpawl_all_NI(trange+1, 3, j, i), '.c', 'MarkerSize', 12);
                        plot(cpawr_all_NI(trange+1, 1, j, i), cpawr_all_NI(trange+1, 3, j, i), '.g', 'MarkerSize', 12);
                    end
                    if ~any(isnan([ctop_all_NI(trange+1, [1 3], j, i) cbottom_all_NI(trange+1, [1 3], j, i)]))
                        plot([ctop_all_NI(trange+1, 1, j, i) cbottom_all_NI(trange+1, 1, j, i)], [ctop_all_NI(trange+1, 3, j, i) cbottom_all_NI(trange+1, 3, j, i)], '-k');
                        plot(ctop_all_NI(trange+1, 1, j, i), ctop_all_NI(trange+1, 3, j, i), '.', 'Color', [1 1 0], 'MarkerSize', 12);
                        plot(cbottom_all_NI(trange+1, 1, j, i), cbottom_all_NI(trange+1, 3, j, i), '.', 'Color', [1 1 0], 'MarkerSize', 12);
                    end
                end
                axis equal;
                xlabel('X (mm)');
                ylabel('Z (mm)');
                set(gca, 'TickDir', 'out', 'FontSize', 12, 'ButtonDownFcn', @extract_figure);
                
                subplot(1, 2, 2);
                hold on;
                for j = 1:nadj_NI(i)
                    if ~any(isnan([cnose_all_NI(trange+1, [1 2], j, i) cpawl_all_NI(trange+1, [1 2], j, i) cpawr_all_NI(trange+1, [1 2], j, i)]))
                        plot([cnose_all_NI(trange+1, 1, j, i) cpawl_all_NI(trange+1, 1, j, i)], [cnose_all_NI(trange+1, 2, j, i) cpawl_all_NI(trange+1, 2, j, i)], '-k');
                        plot([cnose_all_NI(trange+1, 1, j, i) cpawr_all_NI(trange+1, 1, j, i)], [cnose_all_NI(trange+1, 2, j, i) cpawr_all_NI(trange+1, 2, j, i)], '-k');
                        plot([cpawl_all_NI(trange+1, 1, j, i) cpawr_all_NI(trange+1, 1, j, i)], [cpawl_all_NI(trange+1, 2, j, i) cpawr_all_NI(trange+1, 2, j, i)], '-k');
                        plot(cnose_all_NI(trange+1, 1, j, i), cnose_all_NI(trange+1, 2, j, i), '.b', 'MarkerSize', 12);
                        plot(cpawl_all_NI(trange+1, 1, j, i), cpawl_all_NI(trange+1, 2, j, i), '.c', 'MarkerSize', 12);
                        plot(cpawr_all_NI(trange+1, 1, j, i), cpawr_all_NI(trange+1, 2, j, i), '.g', 'MarkerSize', 12);
                    end
                    if ~any(isnan([ctop_all_NI(trange+1, [1 2], j, i) cbottom_all_NI(trange+1, [1 2], j, i)]))
                        plot([ctop_all_NI(trange+1, 1, j, i) cbottom_all_NI(trange+1, 1, j, i)], [ctop_all_NI(trange+1, 2, j, i) cbottom_all_NI(trange+1, 2, j, i)], '-k');
                        plot(ctop_all_NI(trange+1, 1, j, i), ctop_all_NI(trange+1, 2, j, i), '.', 'Color', [1 1 0], 'MarkerSize', 12);
                        plot(cbottom_all_NI(trange+1, 1, j, i), cbottom_all_NI(trange+1, 2, j, i), '.', 'Color', [1 1 0], 'MarkerSize', 12);
                    end
                end
                axis equal;
                xlabel('X (mm)');
                ylabel('Y (mm)');
                set(gca, 'TickDir', 'out', 'FontSize', 12, 'ButtonDownFcn', @extract_figure);
            end
        end
    end
end