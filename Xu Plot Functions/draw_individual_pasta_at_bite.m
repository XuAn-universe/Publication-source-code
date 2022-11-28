function [orientation_NI, orientation_I] = draw_individual_pasta_at_bite(ctop_all_NI, cbottom_all_NI, ctop_all_I, cbottom_all_I, ID_NI, ID_I, trange, i, biteinfo_NI, biteinfo_I)
%*---------------------------------------------------------------------*
% Huang Lab
% Duke University
% Author : Xu An, Dec 2021
% xu.an@duke.edu
% Version: 1.0
%*---------------------------------------------------------------------*
figure;
hold on;
if ~isempty(ID_NI)
    for j = ID_NI
        plot3([ctop_all_NI(trange+1, 1, j, i) cbottom_all_NI(trange+1, 1, j, i)], [ctop_all_NI(trange+1, 2, j, i) cbottom_all_NI(trange+1, 2, j, i)],...
            [ctop_all_NI(trange+1, 3, j, i) cbottom_all_NI(trange+1, 3, j, i)], '-', 'Marker', '.', 'Color', [0.93,0.69,0.13], 'MarkerSize', 12, 'ButtonDownFcn', @displayfilename, 'UserData', biteinfo_NI{j});
    end
end
if ~isempty(ID_I)
    for j = ID_I
        plot3([ctop_all_I(trange+1, 1, j, i) cbottom_all_I(trange+1, 1, j, i)], [ctop_all_I(trange+1, 2, j, i) cbottom_all_I(trange+1, 2, j, i)],...
            [ctop_all_I(trange+1, 3, j, i) cbottom_all_I(trange+1, 3, j, i)], '-', 'Marker', '.', 'Color', [0.92,0.92,0.45], 'MarkerSize', 12, 'ButtonDownFcn', @displayfilename, 'UserData', biteinfo_I{j});
    end
end
axis equal;
xlabel('X (mm)');
ylabel('Y (mm)');
zlabel('Z (mm)');
set(gca, 'TickDir', 'out', 'FontSize', 12, 'ButtonDownFcn', @extract_figure);
view(3);

alpha_NI = [];
beta_NI = [];
gamma_NI = [];
alpha_I = [];
beta_I = [];
gamma_I = [];
anglexy_NI = [];
anglexz_NI = [];
angleyz_NI = [];
anglexy_I = [];
anglexz_I = [];
angleyz_I = [];
if ~isempty(ID_NI)
    dx = ctop_all_NI(trange+1, 1, ID_NI, i)-cbottom_all_NI(trange+1, 1, ID_NI, i);
    dy = ctop_all_NI(trange+1, 2, ID_NI, i)-cbottom_all_NI(trange+1, 2, ID_NI, i);
    dz = ctop_all_NI(trange+1, 3, ID_NI, i)-cbottom_all_NI(trange+1, 3, ID_NI, i);
    dx = squeeze(dx);
    dy = squeeze(dy);
    dz = squeeze(dz);
    alpha_NI = atand(abs(dz./vecnorm([dx dy], 2, 2)));
    beta_NI = atand(abs(dy./vecnorm([dx dz], 2, 2)));
    gamma_NI = atand(abs(dx./vecnorm([dz dy], 2, 2)));
    
    anglexy_NI = atand(dy./dx);
    anglexy_NI(anglexy_NI < 0) = anglexy_NI(anglexy_NI < 0)+180;
    anglexz_NI = atand(dz./dx);
    anglexz_NI(anglexz_NI < 0) = anglexz_NI(anglexz_NI < 0)+180;
    angleyz_NI = atand(dz./dy);
    angleyz_NI(angleyz_NI < 0) = anglexz_NI(angleyz_NI < 0)+180;
end
if ~isempty(ID_I)
    dx = ctop_all_I(trange+1, 1, ID_I, i)-cbottom_all_I(trange+1, 1, ID_I, i);
    dy = ctop_all_I(trange+1, 2, ID_I, i)-cbottom_all_I(trange+1, 2, ID_I, i);
    dz = ctop_all_I(trange+1, 3, ID_I, i)-cbottom_all_I(trange+1, 3, ID_I, i);
    dx = squeeze(dx);
    dy = squeeze(dy);
    dz = squeeze(dz);
    alpha_I = atand(abs(dz./vecnorm([dx dy], 2, 2)));
    beta_I = atand(abs(dy./vecnorm([dx dz], 2, 2)));
    gamma_I = atand(abs(dx./vecnorm([dz dy], 2, 2)));
    
    anglexy_I = atand(dy./dx);
    anglexy_I(anglexy_I < 0) = anglexy_I(anglexy_I < 0)+180;
    anglexz_I = atand(dz./dx);
    anglexz_I(anglexz_I < 0) = anglexz_I(anglexz_I < 0)+180;
    angleyz_I = atand(dz./dy);
    angleyz_I(angleyz_I < 0) = anglexz_I(angleyz_I < 0)+180;
end
orientation_NI.alpha = alpha_NI;
orientation_NI.beta = beta_NI;
orientation_NI.gamma = gamma_NI;
orientation_I.alpha = alpha_I;
orientation_I.beta = beta_I;
orientation_I.gamma = gamma_I;
orientation_NI.anglexy = anglexy_NI;
orientation_NI.anglexz = anglexz_NI;
orientation_NI.angleyz = angleyz_NI;
orientation_I.anglexy = anglexy_I;
orientation_I.anglexz = anglexz_I;
orientation_I.angleyz = angleyz_I;

stp = 2; % default is 2
edges = 0:stp:90;
figure;
subplot(3, 1, 1);
draw_orientations(alpha_NI, alpha_I, edges, 'Orientation XY');

subplot(3, 1, 2);
draw_orientations(beta_NI, beta_I, edges, 'Orientation XZ');

subplot(3, 1, 3);
draw_orientations(gamma_NI, gamma_I, edges, 'Orientation YZ');

edges = 0:stp:180;
figure;
subplot(3, 1, 1);
draw_orientations(anglexy_NI, anglexy_I, edges, 'Orientation XY');

subplot(3, 1, 2);
draw_orientations(anglexz_NI, anglexz_I, edges, 'Orientation XZ');

subplot(3, 1, 3);
draw_orientations(angleyz_NI, angleyz_I, edges, 'Orientation YZ');

    function displayfilename(src, eventdata)
        htext = text(eventdata.IntersectionPoint(1), eventdata.IntersectionPoint(2), src.UserData, 'FontSize', 12, 'Color', 'r');
        pause(2);
        try
            delete(htext)
            clear htext;
        end
    end
end