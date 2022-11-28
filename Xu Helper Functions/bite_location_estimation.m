function [a, b, c, rsquare] = bite_location_estimation(ctop, cbottom, pcolor, method)
%*---------------------------------------------------------------------*
% Huang Lab
% Duke University
% Author : Xu An, Dec 2021
% xu.an@duke.edu
% Version: 1.0
%*---------------------------------------------------------------------*
a = nan;
b = nan;
c = nan;
rsquare = nan;
if size(ctop, 3) <= 1
    return;
end
ctop = squeeze(ctop);
cbottom = squeeze(cbottom);
if ~method
    % assume coordination of bite location is (a, b)
    % (y-b)/(x-a) = k; y-kx = -k*a+b;
    k = (ctop(2, :)-cbottom(2, :))./(ctop(1, :)-cbottom(1, :));
    ynew = ctop(2, :)-k.*ctop(1, :);
    xnew = -k;
    [f, gof] = fit(xnew', ynew', 'poly1', 'Robust', 'Bisquare');
    coefficientValues = coeffvalues(f);
    a = coefficientValues(1);
    b = coefficientValues(2);
    rsquare = gof.rsquare;
    hold on;
    plot(xnew, ynew, 'o', 'Color', pcolor);
    plot([min(xnew) max(xnew)], [a*min(xnew)+b a*max(xnew)+b], '-', 'Color', pcolor);
else
    % k1 = (z2-z1)/(x2-x1); k2 = (y2-y1)/(x2-x1)
    % z-c = k1*(x-a); y-b = k2*(x-a);
    % z-k1*x = -k1*a+c; y-k2*x = -k2*a+b;
    % let y1 = z-k1*x, x1 = -k1, y2 = y-k2*x, x2 = -k2
    % y1 = a*x1+c; y2 = a*x2+b;
    k1 = (ctop(3, :)-cbottom(3, :))./(ctop(1, :)-cbottom(1, :));
    k2 = (ctop(2, :)-cbottom(2, :))./(ctop(1, :)-cbottom(1, :));
    y1 = ctop(3, :)-k1.*ctop(1, :);
    x1 = -k1;
    y2 = ctop(2, :)-k2.*ctop(1, :);
    x2 = -k2;
    y1 = y1';
    x1 = x1';
    y2 = y2';
    x2 = x2';
    xdata = horzcat(x1,x2);
    ydata = horzcat(y1,y2);
    fun = @(b, xdata) horzcat(b(1)*x1+b(3), b(1)*x2+b(2));
    b0 = [4.5; -9; -5.5];
    options = optimoptions('lsqcurvefit', 'MaxFunctionEvaluations', 1000, 'MaxIterations', 1000);
    lb = [];
    ub = [];
    [x, resnorm] = lsqcurvefit(fun, b0, xdata, ydata, lb, ub, options);
    a = x(1);
    b = x(2);
    c = x(3);
    rsquare = 1-resnorm/(sum((y1-mean(y1)).^2)+sum((y2-mean(y2)).^2));
    subplot(1, 2, 1);
    hold on;
    plot(x1, y1, 'o', 'Color', pcolor);
    plot([min(x1) max(x1)], [a*min(x1)+c a*max(x1)+c], '-', 'Color', pcolor);
    title('XZ plane', 'FontSize', 12);
    subplot(1, 2, 2);
    hold on;
    plot(x2, y2, 'o', 'Color', pcolor);
    plot([min(x2) max(x2)], [a*min(x2)+b a*max(x2)+b], '-', 'Color', pcolor);
    title('XY plane', 'FontSize', 12);
end