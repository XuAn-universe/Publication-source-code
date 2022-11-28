function xyz_new = reestimate_ycoordinate_pasta(xyz, yall, zall)
%*---------------------------------------------------------------------*
% Huang Lab
% Duke University
% Author : Xu An, Dec 2021
% xu.an@duke.edu
% Version: 1.0
%*---------------------------------------------------------------------*
xyz_new = xyz;
ynew = [];
for i = [1 4]
    x = zall(i:i+2);
    y = yall(i:i+2);
    temp = isnan(x) | isnan(y);
    x(temp) = [];
    y(temp) = [];
    if numel(x) > 1
        try
            xtemp = ones(numel(x), 2);
            xtemp(:, 1) = x';
            coefficientValues = (xtemp'*xtemp)^(-1)*xtemp'*y';
        catch
            f = fit(x', y', 'poly1');
            coefficientValues = coeffvalues(f);
        end
        a = coefficientValues(1);
        b = coefficientValues(2);
        ynew(end+1) = a*xyz(3)+b;
    end
end
if ~isempty(ynew)
    xyz_new(2) = mean(ynew);
end