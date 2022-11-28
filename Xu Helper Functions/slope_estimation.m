function slope = slope_estimation(x, y)
%*---------------------------------------------------------------------*
% Huang Lab
% Duke University
% Author : Xu An, Dec 2021
% xu.an@duke.edu
% Version: 1.0
%*---------------------------------------------------------------------*
% ncolumn = size(x, 2);
% slope = [];
% for i = 1:ncolumn-1
%     slope = [slope (x(:, i+1:end)-x(:, i))./(y(:, i+1:end)-y(:, i))];
% end
slope = NaN(size(x, 1), 1);
for i = 1:size(x, 1)
    xrow = x(i, :);
    yrow = y(i, :);
    nanid = isnan(xrow) | isnan(yrow);
    xrow(nanid) = [];
    yrow(nanid) = [];
    % slope is y vs x, not x vs y!
    if numel(xrow) > 1
        xrow = xrow';
        yrow = yrow';
        try
            yrow_new = ones(numel(yrow), 2);
            yrow_new(:, 1) = yrow;
            coefficientValues = (yrow_new'*yrow_new)^(-1)*yrow_new'*xrow;
        catch
            f = fit(yrow, xrow, 'poly1');
            coefficientValues = coeffvalues(f);
        end
        slope(i) = coefficientValues(1);
    end
end