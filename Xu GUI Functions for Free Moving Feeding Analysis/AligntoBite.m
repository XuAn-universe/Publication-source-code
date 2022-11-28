function AligntoBite(app, Exp_Path, FrameRate, nmedian)
trial = app.TrialsListBox.Value;
if numel(trial) ~= 1
    helpdlg('You can only choose one trial');
    return;
end

% get bite events
audiolocation = Exp_Path(1:end-7);
try
    temp = load([audiolocation '\Detected_Bite_Events.mat']);
    Bite_events = temp.Audio_analysis;
    bite_timestamps = Bite_events(trial).time_bites;
    bite_amplitudes = Bite_events(trial).amplitude_bites/(max(Bite_events(trial).amplitude_bites));
catch
    bite_timestamps = [];
    bite_amplitudes = [];
end
if isempty(bite_timestamps)
    helpdlg('You must detect the bites first');
    return;
end

if app.sOnTwiceButton.Value
    ninhibition = 2;
    tinhibition = [4 13; 8 17];
elseif app.NoLightButton.Value || app.WholeTrialOnButton.Value || app.BiteDeviceButton.Value
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

if app.sDelay4sOnButton.Value
    bite_timestamps(bite_timestamps <= tinhibition(1) | bite_timestamps > tinhibition(2)) = [];
    if isempty(bite_timestamps)
        return;
    end
end

% process trajectories
if app.TrackingDataCheckBox.Value
    label_table = table2array(app.UITable.Data);
else
    helpdlg('''Tracking Data'' needs to be checked');
    return;
end

[~, ~, ~, ~, ~, xnose, ynose, znose, ~, ~, laserstart, laserstop] = trajectory_postprocessing(3, Exp_Path, trial, label_table, nmedian, FrameRate);
[~, ~, ~, ~, ~, xpawl, ypawl, zpawl, ~, ~, ~, ~] = trajectory_postprocessing(10, Exp_Path, trial, label_table, nmedian, FrameRate);
[~, ~, ~, ~, ~, xpawr, ypawr, zpawr, ~, ~, ~, ~] = trajectory_postprocessing(16, Exp_Path, trial, label_table, nmedian, FrameRate);
[y_topPG1, y_topPG3, z_topPG1, ~, z_topPG3, xtop, ytop, ztop, ~, ~, ~, ~] = trajectory_postprocessing(21, Exp_Path, trial, label_table, nmedian, FrameRate);
[y_bottomPG1, y_bottomPG3, z_bottomPG1, ~, z_bottomPG3, xbottom, ybottom, zbottom, ~, ~, ~, ~] = trajectory_postprocessing(22, Exp_Path, trial, label_table, nmedian, FrameRate);
[y_centerPG1, y_centerPG3, z_centerPG1, ~, z_centerPG3, xcenter, ycenter, zcenter, ~, ~, ~, ~] = trajectory_postprocessing(33, Exp_Path, trial, label_table, nmedian, FrameRate);

if analyze_method == 3
    [~, ~, ~, ~, ~, xanklel, yanklel, zanklel, ~, ~, ~, ~] = trajectory_postprocessing(17, Exp_Path, trial, label_table, nmedian, FrameRate);
    [~, ~, ~, ~, ~, xankler, yankler, zankler, ~, ~, ~, ~] = trajectory_postprocessing(18, Exp_Path, trial, label_table, nmedian, FrameRate);
    corigin = [(xanklel(:, 1)+xankler(:, 1))/2 (yanklel(:, 1)+yankler(:, 1))/2 (zanklel(:, 1)+zankler(:, 1))/2];
    vankleRL = [xanklel(:, 1)-xankler(:, 1) yanklel(:, 1)-yankler(:, 1)];
    theta = atan(vankleRL(:, 2)./vankleRL(:, 1));
    theta(vankleRL(:, 1) < 0) = pi+theta(vankleRL(:, 1) < 0);
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
    ctop = [xtop(:, 1) ytop(:, 1) ztop(:, 1)];
    cbottom = [xbottom(:, 1) ybottom(:, 1) zbottom(:, 1)];
    ccenter = [xcenter(:, 1) ycenter(:, 1) zcenter(:, 1)];
    for i = 1:size(xnose, 1)
        % re-estimate pasta y coordinate based on z coordinate
        y_pasta = [y_topPG1(i) y_bottomPG1(i) y_centerPG1(i) y_topPG3(i) y_bottomPG3(i) y_centerPG3(i)];
        z_pasta = [z_topPG1(i) z_bottomPG1(i) z_centerPG1(i) z_topPG3(i) z_bottomPG3(i) z_centerPG3(i)];
        ctop(i, :) = reestimate_ycoordinate_pasta(ctop(i, :), y_pasta, z_pasta);
        cbottom(i, :) = reestimate_ycoordinate_pasta(cbottom(i, :), y_pasta, z_pasta);
        ccenter(i, :) = reestimate_ycoordinate_pasta(ccenter(i, :), y_pasta, z_pasta);
    end
    ctopT = [ctop(:, 1)-corigin(:, 1) ctop(:, 2)-corigin(:, 2) ctop(:, 3)-corigin(:, 3)];
    ctop(:, 1) = ctopT(:, 1).*cos(theta)+ctopT(:, 2).*sin(theta);
    ctop(:, 2) = ctopT(:, 2).*cos(theta)-ctopT(:, 1).*sin(theta);
    ctop(:, 3) = ctopT(:, 3);
    cbottomT = [cbottom(:, 1)-corigin(:, 1) cbottom(:, 2)-corigin(:, 2) cbottom(:, 3)-corigin(:, 3)];
    cbottom(:, 1) = cbottomT(:, 1).*cos(theta)+cbottomT(:, 2).*sin(theta);
    cbottom(:, 2) = cbottomT(:, 2).*cos(theta)-cbottomT(:, 1).*sin(theta);
    cbottom(:, 3) = cbottomT(:, 3);
    ccenterT = [ccenter(:, 1)-corigin(:, 1) ccenter(:, 2)-corigin(:, 2) ccenter(:, 3)-corigin(:, 3)];
    ccenter(:, 1) = ccenterT(:, 1).*cos(theta)+ccenterT(:, 2).*sin(theta);
    ccenter(:, 2) = ccenterT(:, 2).*cos(theta)-ccenterT(:, 1).*sin(theta);
    ccenter(:, 3) = ccenterT(:, 3);
    
    t = (1:size(xnose, 1))'/FrameRate;
    
    lc1 = ones(trange+1, 3);
    lc1(:, 1) = 0.5;
    lc1(:, 1) = linspace(0, 0.85, trange+1)';
    lc2 = ones(trange+1, 3);
    lc2(:, 1) = 0.3;
    lc2(:, 1) = linspace(0, 0.85, trange+1)';
    lc3 = ones(trange+1, 3);
    lc3(:, 1) = 0;
    lc3(:, 1) = linspace(0, 0.85, trange+1)';
    lc4 = ones(trange+1, 3);
    lc4(:, 1) = 0.12;
    lc4(:, 1) = linspace(0, 0.85, trange+1)';
    
    for i = 1:numel(bite_timestamps)
        timestamp = bite_timestamps(i);
        [~, id] = min(abs(t-timestamp));
        if id-trange >= 1
            figure;
            subplot(1, 2, 1);
            hold on;
            for j = id-trange:id
                if ~any(isnan([cnose(j, [1 3]) cpawl(j, [1 3]) cpawr(j, [1 3]) ctop(j, [1 3]) cbottom(j, [1 3])]))
                    plot([cnose(j, 1) cpawl(j, 1)], [cnose(j, 3) cpawl(j, 3)], '-', 'Color', hsv2rgb(lc1(j-id+trange+1, :)));
                    plot([cnose(j, 1) cpawr(j, 1)], [cnose(j, 3) cpawr(j, 3)], '-', 'Color', hsv2rgb(lc2(j-id+trange+1, :)));
                    plot([cpawr(j, 1) cpawl(j, 1)], [cpawr(j, 3) cpawl(j, 3)], '-', 'Color', hsv2rgb(lc3(j-id+trange+1, :)));
                    plot([ctop(j, 1) cbottom(j, 1)], [ctop(j, 3) cbottom(j, 3)], '-', 'Color', hsv2rgb(lc4(j-id+trange+1, :)));
                    plot(cnose(j, 1), cnose(j, 3), '.b', 'MarkerSize', 12);
                    plot(cpawl(j, 1), cpawl(j, 3), '.c', 'MarkerSize', 12);
                    plot(cpawr(j, 1), cpawr(j, 3), '.g', 'MarkerSize', 12);
                    plot(ctop(j, 1), ctop(j, 3), '.k', 'MarkerSize', 12);
                    plot(cbottom(j, 1), cbottom(j, 3), '.k', 'MarkerSize', 12);
                    %             plot(ccenter(j, 1), ccenter(j, 3), 'ko');
                end
            end
            axis equal;
            xlabel('X (mm)');
            ylabel('Z (mm)');
            set(gca, 'TickDir', 'out', 'FontSize', 12, 'ButtonDownFcn', @extract_figure);
            
            subplot(1, 2, 2);
            hold on;
            for j = id-trange:id
                if ~any(isnan([cnose(j, [1 2]) cpawl(j, [1 2]) cpawr(j, [1 2]) ctop(j, [1 2]) cbottom(j, [1 2])]))
                    plot([cnose(j, 1) cpawl(j, 1)], [cnose(j, 2) cpawl(j, 2)], '-', 'Color', hsv2rgb(lc1(j-id+trange+1, :)));
                    plot([cnose(j, 1) cpawr(j, 1)], [cnose(j, 2) cpawr(j, 2)], '-', 'Color', hsv2rgb(lc2(j-id+trange+1, :)));
                    plot([cpawr(j, 1) cpawl(j, 1)], [cpawr(j, 2) cpawl(j, 2)], '-', 'Color', hsv2rgb(lc3(j-id+trange+1, :)));
                    plot([ctop(j, 1) cbottom(j, 1)], [ctop(j, 2) cbottom(j, 2)], '-', 'Color', hsv2rgb(lc4(j-id+trange+1, :)));
                    plot(cnose(j, 1), cnose(j, 2), '.b', 'MarkerSize', 12);
                    plot(cpawl(j, 1), cpawl(j, 2), '.c', 'MarkerSize', 12);
                    plot(cpawr(j, 1), cpawr(j, 2), '.g', 'MarkerSize', 12);
                    plot(ctop(j, 1), ctop(j, 2), '.k', 'MarkerSize', 12);
                    plot(cbottom(j, 1), cbottom(j, 2), '.k', 'MarkerSize', 12);
                    %             plot(ccenter(j, 1), ccenter(j, 2), 'ko');
                end
            end
            axis equal;
            xlabel('X (mm)');
            ylabel('Y (mm)');
            set(gca, 'TickDir', 'out', 'FontSize', 12, 'ButtonDownFcn', @extract_figure);
        end
    end
elseif analyze_method == 2
    t = (1:size(xnose, 1))'/FrameRate;
    znose_all = nan(trange*2+1, numel(bite_timestamps));
    zpawl_all = nan(size(znose_all));
    zpawr_all = nan(size(znose_all));
    ts = (-trange:1:trange)'/FrameRate;
    time_range = [ts(1) ts(end)];
    for i = 1:numel(bite_timestamps)
        timestamp = bite_timestamps(i);
        [~, id] = min(abs(t-timestamp));
        if id-trange >= 1
            znose_all(:, i) = znose(id-trange:id+trange, 1);
            zpawl_all(:, i) = zpawl(id-trange:id+trange, 1);
            zpawr_all(:, i) = zpawr(id-trange:id+trange, 1);
            
            figure;
            hold on;
            plot(ts, znose(id-trange:id+trange, 1), '-b');
            plot(ts, zpawl(id-trange:id+trange, 1), '-c');
            plot(ts, zpawr(id-trange:id+trange, 1), '-g');
        end
    end
    
    figure;
    plot_tj_multitrial(time_range, ts, znose_all, [], [], [], [], '', '');
    plot_tj_multitrial(time_range, ts, zpawl_all, [], [], [], [], '', '');
    plot_tj_multitrial(time_range, ts, zpawr_all, [], [], [], [], 'Z (mm)', '');
    ylim([min([znose_all(:); zpawl_all(:); zpawr_all(:)]) max([znose_all(:); zpawl_all(:); zpawr_all(:)])]);
elseif analyze_method == 1
    [~, ~, ~, ~, ~, xeyel, yeyel, zeyel, ~, ~, ~, ~] = trajectory_postprocessing(1, Exp_Path, trial, label_table, nmedian, FrameRate);
    [~, ~, ~, ~, ~, xeyer, yeyer, zeyer, ~, ~, ~, ~] = trajectory_postprocessing(2, Exp_Path, trial, label_table, nmedian, FrameRate);
%     corigin = [xnose(:, 1) ynose(:, 1) znose(:, 1)];
%     % translation, new coordinate origin is nose
    corigin = [xeyer(:, 1) yeyer(:, 1) zeyer(:, 1)];
    % translation, new coordinate origin is right eye
    ceyelT = [xeyel(:, 1)-corigin(:, 1) yeyel(:, 1)-corigin(:, 2) zeyel(:, 1)-corigin(:, 3)];
%     ceyerT = [xeyer(:, 1)-corigin(:, 1) yeyer(:, 1)-corigin(:, 2) zeyer(:, 1)-corigin(:, 3)];
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
    
    ceyel_all = nan(3, size(xnose, 1));
%     ceyer_all = nan(3, size(xnose, 1));
    cnose_all = nan(3, size(xnose, 1));
    cpawl_all = nan(3, size(xnose, 1));
    cpawr_all = nan(3, size(xnose, 1));
    ctop_all = nan(3, size(xnose, 1));
    cbottom_all = nan(3, size(xnose, 1));
    ccenter_all = nan(3, size(xnose, 1));
    for i = 1:size(xnose, 1)
        ceyel = ceyelT(i, :);
%         ceyer = ceyerT(i, :);
        cnose = cnoseT(i, :);
%         if ~any(isnan([ceyel ceyer]))
        if ~any(isnan([ceyel cnose]))
%             a = norm(ceyer);
            a = norm(cnose);
%             b = norm(ceyel);
            b = norm(ceyel-cnose);
%             c = norm(ceyel-ceyer);
            c = norm(ceyel);
%             y = sqrt(a^2*c^2-dot(-ceyer, ceyel-ceyer)^2)/c; % new y coordinate of left & right eyes
            y = -sqrt(a^2*c^2-dot(cnose, ceyel)^2)/c; % new y coordinate of nose
%             x1 = -dot(-ceyer, ceyel-ceyer)/c; % new x coordinate of right eye
            x1 = dot(cnose, ceyel)/c; % new x coordinate of nose
%             x2 = dot(ceyel, ceyel-ceyer)/c; % new x coordinate of left eye
            x2 = c; % new x coordinate of left eye
            if 0
                if x1 < -10 || x2 < 2.5 || x1 > -2.5 || x2 > 10 || x2-x1 < 5 || y < 11.5 || y > 13
                    continue;
                end
            end
            % standard othogonal bases of the new coordinate system
%             u = (ceyel-ceyer)/c;
%             v = (ceyer+(ceyel-ceyer)*(-x1)/c)/y;
%             w = cross(ceyel, ceyer)/norm(cross(ceyel, ceyer));
            u = ceyel/c;
            v = (-cnose+ceyel/c*x1)/(-y);
            w = cross(cnose, ceyel)/norm(cross(cnose, ceyel));
            
            cpawl = cpawlT(i, :);
            cpawr = cpawrT(i, :);
            ctop = ctopT(i, :);
            cbottom = cbottomT(i, :);
            ccenter = ccenterT(i, :);
            % re-estimate pasta y coordinate based on z coordinate
            y_pasta = y_pastaT(i, :);
            z_pasta = z_pastaT(i, :);
            ctop = reestimate_ycoordinate_pasta(ctop, y_pasta, z_pasta);
            cbottom = reestimate_ycoordinate_pasta(cbottom, y_pasta, z_pasta);
            ccenter = reestimate_ycoordinate_pasta(ccenter, y_pasta, z_pasta);
            
            % compute the new coordinates
            ceyel_all(:, i) = [u; v; w]*ceyel';
%             ceyer_all(:, i) = [u; v; w]*ceyer';
            cnose_all(:, i) = [u; v; w]*cnose';
            cpawl_all(:, i) = [u; v; w]*cpawl';
            cpawr_all(:, i) = [u; v; w]*cpawr';
            ctop_all(:, i) = [u; v; w]*ctop';
            cbottom_all(:, i) = [u; v; w]*cbottom';
            ccenter_all(:, i) = [u; v; w]*ccenter';
        end
    end
    
    t = (1:size(xnose, 1))'/FrameRate;
    
    lc1 = ones(trange+1, 3);
    lc1(:, 1) = linspace(0, 0.85, trange+1)';
    lc2 = ones(trange+1, 3);
    lc2(:, 1) = linspace(0, 0.85, trange+1)';
    lc3 = ones(trange+1, 3);
    lc3(:, 1) = linspace(0, 0.85, trange+1)';
    lc4 = ones(trange+1, 3);
    lc4(:, 1) = linspace(0, 0.85, trange+1)';
    lc5 = ones(trange+1, 3);
    lc5(:, 1) = linspace(0, 0.85, trange+1)';
    lc6 = ones(trange+1, 3);
    lc6(:, 1) = linspace(0, 0.85, trange+1)';
    lc7 = ones(trange+1, 3);
    lc7(:, 1) = linspace(0, 0.85, trange+1)';
    
    for i = 1:numel(bite_timestamps)
        timestamp = bite_timestamps(i);
        [~, id] = min(abs(t-timestamp));
        if id-trange >= 1
            figure;
            subplot(1, 2, 1);
            hold on;
            for j = id-trange:id
%                 plot([ceyer_all(1, j) ceyel_all(1, j)], [ceyer_all(3, j) ceyel_all(3, j)], '-', 'Color', hsv2rgb(lc1(j-id+trange+1, :)));
                plot([0 ceyel_all(1, j)], [0 ceyel_all(3, j)], '-', 'Color', hsv2rgb(lc1(j-id+trange+1, :)));
%                 plot([0 ceyel_all(1, j)], [0 ceyel_all(3, j)], '-', 'Color', hsv2rgb(lc2(j-id+trange+1, :)));
                plot([cnose_all(1, j) ceyel_all(1, j)], [cnose_all(3, j) ceyel_all(3, j)], '-', 'Color', hsv2rgb(lc2(j-id+trange+1, :)));
%                 plot([0 ceyer_all(1, j)], [0 ceyer_all(3, j)], '-', 'Color', hsv2rgb(lc3(j-id+trange+1, :)));
                plot([cnose_all(1, j) 0], [cnose_all(3, j) 0], '-', 'Color', hsv2rgb(lc3(j-id+trange+1, :)));
                plot(ceyel_all(1, j), ceyel_all(3, j), '.', 'Color', [0 0 0], 'MarkerSize', 12);
%                 plot(ceyer_all(1, j), ceyer_all(3, j), '.', 'Color', [0 0 0], 'MarkerSize', 12);
                plot(0, 0, '.', 'Color', [0 0 0], 'MarkerSize', 12);
                
                if ~any(isnan([cpawl_all(1, j) cpawl_all(3, j) cpawr_all(1, j) cpawr_all(3, j)]))
%                     plot([0 cpawl_all(1, j)], [0 cpawl_all(3, j)], '-', 'Color', hsv2rgb(lc4(j-id+trange+1, :)));
%                     plot([0 cpawr_all(1, j)], [0 cpawr_all(3, j)], '-', 'Color', hsv2rgb(lc5(j-id+trange+1, :)));
                    plot([cnose_all(1, j) cpawl_all(1, j)], [cnose_all(3, j) cpawl_all(3, j)], '-', 'Color', hsv2rgb(lc4(j-id+trange+1, :)));
                    plot([cnose_all(1, j) cpawr_all(1, j)], [cnose_all(3, j) cpawr_all(3, j)], '-', 'Color', hsv2rgb(lc5(j-id+trange+1, :)));
                    plot([cpawl_all(1, j) cpawr_all(1, j)], [cpawl_all(3, j) cpawr_all(3, j)], '-', 'Color', hsv2rgb(lc6(j-id+trange+1, :)));
                    
                    plot(cpawl_all(1, j), cpawl_all(3, j), '.', 'Color', [0.00,0.45,0.74], 'MarkerSize', 12);
                    plot(cpawr_all(1, j), cpawr_all(3, j), '.', 'Color', [0.47,0.67,0.19], 'MarkerSize', 12);
                end
                
                if ~any(isnan([ctop_all(1, j) ctop_all(3, j) cbottom_all(1, j) cbottom_all(3, j)]))
                    plot([ctop_all(1, j) cbottom_all(1, j)], [ctop_all(3, j) cbottom_all(3, j)], '-', 'Color', hsv2rgb(lc7(j-id+trange+1, :)));
                    
                    plot(ctop_all(1, j), ctop_all(3, j), '.', 'Color', [0.93,0.69,0.13], 'MarkerSize', 12);
                    plot(cbottom_all(1, j), cbottom_all(3, j), '.', 'Color', [0.93,0.69,0.13], 'MarkerSize', 12);
                end
            end
            axis equal;
            xlabel('X (mm)');
            ylabel('Z (mm)');
            set(gca, 'TickDir', 'out', 'FontSize', 12, 'ButtonDownFcn', @extract_figure);
            
            subplot(1, 2, 2);
            hold on;
            for j = id-trange:id
%                 plot([ceyer_all(1, j) ceyel_all(1, j)], [ceyer_all(2, j) ceyel_all(2, j)], '-', 'Color', hsv2rgb(lc1(j-id+trange+1, :)));
%                 plot([0 ceyel_all(1, j)], [0 ceyel_all(2, j)], '-', 'Color', hsv2rgb(lc2(j-id+trange+1, :)));
%                 plot([0 ceyer_all(1, j)], [0 ceyer_all(2, j)], '-', 'Color', hsv2rgb(lc3(j-id+trange+1, :)));
                plot([0 ceyel_all(1, j)], [0 ceyel_all(2, j)], '-', 'Color', hsv2rgb(lc1(j-id+trange+1, :)));
                plot([cnose_all(1, j) ceyel_all(1, j)], [cnose_all(2, j) ceyel_all(2, j)], '-', 'Color', hsv2rgb(lc2(j-id+trange+1, :)));
                plot([cnose_all(1, j) 0], [cnose_all(2, j) 0], '-', 'Color', hsv2rgb(lc3(j-id+trange+1, :)));
                plot(ceyel_all(1, j), ceyel_all(2, j), '.', 'Color', [0 0 0], 'MarkerSize', 12);
%                 plot(ceyer_all(1, j), ceyer_all(2, j), '.', 'Color', [0 0 0], 'MarkerSize', 12);
                plot(0, 0, '.', 'Color', [0 0 0], 'MarkerSize', 12);
                
                if ~any(isnan([cpawl_all(1, j) cpawl_all(2, j) cpawr_all(1, j) cpawr_all(2, j)]))
%                     plot([0 cpawl_all(1, j)], [0 cpawl_all(2, j)], '-', 'Color', hsv2rgb(lc4(j-id+trange+1, :)));
%                     plot([0 cpawr_all(1, j)], [0 cpawr_all(2, j)], '-', 'Color', hsv2rgb(lc5(j-id+trange+1, :)));
                    plot([cnose_all(1, j) cpawl_all(1, j)], [cnose_all(2, j) cpawl_all(2, j)], '-', 'Color', hsv2rgb(lc4(j-id+trange+1, :)));
                    plot([cnose_all(1, j) cpawr_all(1, j)], [cnose_all(2, j) cpawr_all(2, j)], '-', 'Color', hsv2rgb(lc5(j-id+trange+1, :)));
                    plot([cpawl_all(1, j) cpawr_all(1, j)], [cpawl_all(2, j) cpawr_all(2, j)], '-', 'Color', hsv2rgb(lc6(j-id+trange+1, :)));
                    
                    plot(cpawl_all(1, j), cpawl_all(2, j), '.', 'Color', [0.00,0.45,0.74], 'MarkerSize', 12);
                    plot(cpawr_all(1, j), cpawr_all(2, j), '.', 'Color', [0.47,0.67,0.19], 'MarkerSize', 12);
                end
                
                if ~any(isnan([ctop_all(1, j) ctop_all(2, j) cbottom_all(1, j) cbottom_all(2, j)]))
                    plot([ctop_all(1, j) cbottom_all(1, j)], [ctop_all(2, j) cbottom_all(2, j)], '-', 'Color', hsv2rgb(lc7(j-id+trange+1, :)));
                    
                    plot(ctop_all(1, j), ctop_all(2, j), '.', 'Color', [0.93,0.69,0.13], 'MarkerSize', 12);
                    plot(cbottom_all(1, j), cbottom_all(2, j), '.', 'Color', [0.93,0.69,0.13], 'MarkerSize', 12);
                end
            end
            axis equal;
            xlabel('X (mm)');
            ylabel('Y (mm)');
            set(gca, 'TickDir', 'out', 'FontSize', 12, 'ButtonDownFcn', @extract_figure);
        end
    end
end