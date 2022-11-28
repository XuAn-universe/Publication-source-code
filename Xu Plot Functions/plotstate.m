function plotstate(t, statepath, yrange, pcolor)
%*---------------------------------------------------------------------*
% Huang Lab
% Duke University
% Author : Xu An, Jan 2022
% xu.an@duke.edu
% Version: 1.0
%*---------------------------------------------------------------------*
states = unique(statepath);
hold on;
for i = 1:numel(states)
    temp = diff([0 statepath == states(i) 0]);
    xstart = t(temp == 1);
    xstop = t(find(temp == -1)-1);
    
    patch([xstart'; xstart'; xstop'; xstop'], repmat([yrange(1); yrange(2); yrange(2); yrange(1)], 1, numel(xstart)), pcolor(i, :), 'EdgeColor', 'none', 'FaceAlpha', 0.2);
end
hold off;