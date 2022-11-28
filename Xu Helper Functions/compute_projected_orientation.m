function [orientation_NI, orientation_I] = compute_projected_orientation(ctopNI, cbottomNI, ctopI, cbottomI, pcheck, stp)
%*---------------------------------------------------------------------*
% Huang Lab
% Duke University
% Author : Xu An, Dec 2021
% xu.an@duke.edu
% Version: 1.0
%*---------------------------------------------------------------------*
if nargin < 6
    stp = 10;
end
oxy = atand((ctopNI(2, :)-cbottomNI(2, :))./(ctopNI(1, :)-cbottomNI(1, :)));
oxy(oxy < 0) = oxy(oxy < 0)+180;
% oxy(oxy < 0) = abs(oxy(oxy < 0));
oxz = atand((ctopNI(3, :)-cbottomNI(3, :))./(ctopNI(1, :)-cbottomNI(1, :)));
oxz(oxz < 0) = oxz(oxz < 0)+180;
% oxz(oxz < 0) = abs(oxz(oxz < 0));
oyz = atand((ctopNI(3, :)-cbottomNI(3, :))./(ctopNI(2, :)-cbottomNI(2, :)));
oyz(oyz < 0) = oyz(oyz < 0)+180;
% oyz(oyz < 0) = abs(oyz(oyz < 0));
orientation_NI.xy = oxy;
orientation_NI.xz = oxz;
orientation_NI.yz = oyz;

if ~isempty(ctopI)
    oxy = atand((ctopI(2, :)-cbottomI(2, :))./(ctopI(1, :)-cbottomI(1, :)));
    oxy(oxy < 0) = oxy(oxy < 0)+180;
%     oxy(oxy < 0) = abs(oxy(oxy < 0));
    oxz = atand((ctopI(3, :)-cbottomI(3, :))./(ctopI(1, :)-cbottomI(1, :)));
    oxz(oxz < 0) = oxz(oxz < 0)+180;
%     oxz(oxz < 0) = abs(oxz(oxz < 0));
    oyz = atand((ctopI(3, :)-cbottomI(3, :))./(ctopI(2, :)-cbottomI(2, :)));
    oyz(oyz < 0) = oyz(oyz < 0)+180;
%     oyz(oyz < 0) = abs(oyz(oyz < 0));
    orientation_I.xy = oxy;
    orientation_I.xz = oxz;
    orientation_I.yz = oyz;
else
    orientation_I.xy = [];
    orientation_I.xz = [];
    orientation_I.yz = [];
end

if pcheck
    edges = 0:stp:180;
    figure;
    subplot(3, 1, 1);
    draw_orientations(orientation_NI.xz, orientation_I.xz, edges, 'Orientation XZ');
    subplot(3, 1, 2);
    draw_orientations(orientation_NI.xy, orientation_I.xy, edges, 'Orientation XY');
    subplot(3, 1, 3);
    draw_orientations(orientation_NI.yz, orientation_I.yz, edges, 'Orientation YZ');
end
end