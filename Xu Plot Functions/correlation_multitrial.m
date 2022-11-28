function correlation_multitrial(var1_NI, var2_NI, var1_I, var2_I, xlabel_text, ylabel_text, title_text, axis_setting)
ID_NI = ~isnan(var1_NI) & ~isnan(var2_NI);
mdl_NI = fitlm(var1_NI(ID_NI), var2_NI(ID_NI))
[rho_NI, pval_NI] = corr(var1_NI(ID_NI), var2_NI(ID_NI))
figure;
plot(var1_NI(ID_NI), var2_NI(ID_NI), 'ok');
hold on;

if ~isempty(var1_I) && ~isempty(var2_I)
    ID_I = ~isnan(var1_I) & ~isnan(var2_I);
    mdl_I = fitlm(var1_I(ID_I), var2_I(ID_I))
    [rho_I, pval_I] = corr(var1_I(ID_I), var2_I(ID_I))
    plot(var1_I(ID_I), var2_I(ID_I), 'og');
end

xnew = sort(var1_NI(ID_NI));
[ypred_NI, yci_NI] = predict(mdl_NI, xnew);
plot(xnew, ypred_NI, '-k');
plot(xnew, yci_NI, '--k');

if ~isempty(var1_I) && ~isempty(var2_I)
    xnew = sort(var1_I(ID_I));
    [ypred_I, yci_I] = predict(mdl_I, xnew);
    plot(xnew, ypred_I, '-', 'Color', [0 0.75 0]);
    plot(xnew, yci_I, '--', 'Color', [0 0.75 0]);
end
if ~isempty(axis_setting)
    axis(axis_setting);
end
grid on;
xlabel(xlabel_text);
ylabel(ylabel_text);
title(title_text);
set(gca, 'FontSize', 12);