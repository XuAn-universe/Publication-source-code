function result = correlation_plot(xcontrol, ycontrol, nleft_control, xinhibition, yinhibition, nleft_inhibition, xlabel_text, ylabel_text, orientation_check)
%*---------------------------------------------------------------------*
% Huang Lab
% Duke University
% Author : Xu An, Jan 2022
% xu.an@duke.edu
% Version: 1.0
%*---------------------------------------------------------------------*
thr = 15;
mksize = 8;
xfit_NI = [];
yfit_NI = [];
xfit_I = [];
yfit_I = [];
if ~isempty(xcontrol) && ~isempty(ycontrol)
    xfit_NI = [xfit_NI; xcontrol(:)];
    yfit_NI = [yfit_NI; ycontrol(:)];
end
if ~isempty(xinhibition) && ~isempty(yinhibition)
    xfit_I = [xfit_I; xinhibition(:)];
    yfit_I = [yfit_I; yinhibition(:)];
end
IDexc = find(isnan(xfit_NI)|isnan(yfit_NI));
nleft_control = nleft_control-sum(IDexc <= nleft_control);
xfit_NI(IDexc) = [];
yfit_NI(IDexc) = [];
IDexc = find(isnan(xfit_I)|isnan(yfit_I));
nleft_inhibition = nleft_inhibition-sum(IDexc <= nleft_inhibition);
xfit_I(IDexc) = [];
yfit_I(IDexc) = [];

rsquare = nan(1, 4); % R_NI, L_NI, R_I, L_I
slope = nan(1, 4); % R_NI, L_NI, R_I, L_I
rho = nan(1, 4); % R_NI, L_NI, R_I, L_I
p = nan(1, 4); % R_NI, L_NI, R_I, L_I

hold on; % left light, right dark
try
    x = xfit_NI(nleft_control+1:end);
    y = yfit_NI(nleft_control+1:end);
    if orientation_check
        y(y > 90) = y(y > 90)-180;
    end
    if numel(x) > thr
        scatter(x, y, mksize, [0 0 0], 'filled');
        [f, gof] = fit(x, y, 'poly1');
        coefficientValues = coeffvalues(f);
        a = coefficientValues(1);
        b = coefficientValues(2);
        plot([min(x) max(x)], [a*min(x)+b a*max(x)+b], '-', 'Color', [0 0 0]);
        
        rsquare(1) = gof.rsquare;
        slope(1) = a;
        [rho(1), p(1)] = corr(x, y);
    end
end
try
    x = xfit_NI(1:nleft_control);
    y = yfit_NI(1:nleft_control);
    if orientation_check
        y(y < 90) = y(y < 90)+180;
    end
    if numel(x) > thr
        scatter(x, y, mksize, [0.8 0.8 0.8], 'filled');
        [f, gof] = fit(x, y, 'poly1');
        coefficientValues = coeffvalues(f);
        a = coefficientValues(1);
        b = coefficientValues(2);
        plot([min(x) max(x)], [a*min(x)+b a*max(x)+b], '-', 'Color', [0.8 0.8 0.8]);
        
        rsquare(2) = gof.rsquare;
        slope(2) = a;
        [rho(2), p(2)] = corr(x, y);
    end
end
try
    x = xfit_I(nleft_inhibition+1:end);
    y = yfit_I(nleft_inhibition+1:end);
    if orientation_check
        y(y > 90) = y(y > 90)-180;
    end
    if numel(x) > thr
        scatter(x, y, mksize, [0 0 1], 'filled');
        [f, gof] = fit(x, y, 'poly1');
        coefficientValues = coeffvalues(f);
        a = coefficientValues(1);
        b = coefficientValues(2);
        plot([min(x) max(x)], [a*min(x)+b a*max(x)+b], '-', 'Color', [0 0 1]);
        
        rsquare(3) = gof.rsquare;
        slope(3) = a;
        [rho(3), p(3)] = corr(x, y);
    end
end
try
    x = xfit_I(1:nleft_inhibition);
    y = yfit_I(1:nleft_inhibition);
    if orientation_check
        y(y < 90) = y(y < 90)+180;
    end
    if numel(x) > thr
        scatter(x, y, mksize, [0.3010 0.7450 0.9330], 'filled');
        [f, gof] = fit(x, y, 'poly1');
        coefficientValues = coeffvalues(f);
        a = coefficientValues(1);
        b = coefficientValues(2);
        plot([min(x) max(x)], [a*min(x)+b a*max(x)+b], '-', 'Color', [0.3010 0.7450 0.9330]);
        
        rsquare(4) = gof.rsquare;
        slope(4) = a;
        [rho(4), p(4)] = corr(x, y);
    end
end
set(gca, 'TickDir', 'in', 'FontSize', 12, 'ButtonDownFcn', @extract_figure);
box off;
xlabel(xlabel_text);
ylabel(ylabel_text);
title(['R_NI=' num2str(rsquare(1), '%.2f') ';L_NI=' num2str(rsquare(2), '%.2f') ';R_I=' num2str(rsquare(3), '%.2f') ';L_I=' num2str(rsquare(4), '%.2f')]);

result.rsquare = rsquare;
result.slope = slope;
result.rho = rho;
result.p = p;