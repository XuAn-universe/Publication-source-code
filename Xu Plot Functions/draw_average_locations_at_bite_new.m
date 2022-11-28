function draw_average_locations_at_bite_new(ceyel_all_NI, cnose_all_NI, cpawl_all_NI, cpawr_all_NI, ctop_all_NI, cbottom_all_NI,...
    ceyel_all_I, cnose_all_I, cpawl_all_I, cpawr_all_I, ctop_all_I, cbottom_all_I, ID_front_NI, ID_top_NI, ID_front_I, ID_top_I, trange, i)
%*---------------------------------------------------------------------*
% Huang Lab
% Duke University
% Author : Xu An, Dec 2021
% xu.an@duke.edu
% Version: 1.0
%*---------------------------------------------------------------------*
figure;
subplot(1, 2, 1);
hold on;
if ~isempty(ID_front_NI)
    plotline(0, ceyel_all_NI(trange+1, 1, ID_front_NI, i), 0, ceyel_all_NI(trange+1, 3, ID_front_NI, i), [0 0 0]);
    plotline(cnose_all_NI(trange+1, 1, ID_front_NI, i), ceyel_all_NI(trange+1, 1, ID_front_NI, i), cnose_all_NI(trange+1, 3, ID_front_NI, i), ceyel_all_NI(trange+1, 3, ID_front_NI, i), [0 0 0]);
    plotline(cnose_all_NI(trange+1, 1, ID_front_NI, i), 0, cnose_all_NI(trange+1, 3, ID_front_NI, i), 0, [0 0 0]);
    
    plotline(cnose_all_NI(trange+1, 1, ID_front_NI, i), cpawl_all_NI(trange+1, 1, ID_front_NI, i), cnose_all_NI(trange+1, 3, ID_front_NI, i), cpawl_all_NI(trange+1, 3, ID_front_NI, i), [0 0 0]);
    plotline(cnose_all_NI(trange+1, 1, ID_front_NI, i), cpawr_all_NI(trange+1, 1, ID_front_NI, i), cnose_all_NI(trange+1, 3, ID_front_NI, i), cpawr_all_NI(trange+1, 3, ID_front_NI, i), [0 0 0]);
    plotline(cpawl_all_NI(trange+1, 1, ID_front_NI, i), cpawr_all_NI(trange+1, 1, ID_front_NI, i), cpawl_all_NI(trange+1, 3, ID_front_NI, i), cpawr_all_NI(trange+1, 3, ID_front_NI, i), [0 0 0]);
    
    plotline(ctop_all_NI(trange+1, 1, ID_front_NI, i), cbottom_all_NI(trange+1, 1, ID_front_NI, i),...
        ctop_all_NI(trange+1, 3, ID_front_NI, i), cbottom_all_NI(trange+1, 3, ID_front_NI, i), [0.93,0.69,0.13]);

    ploterrorbar(ceyel_all_NI(trange+1, 1, ID_front_NI, i), 0, [0 0 0]);
    ploterrorbar(cnose_all_NI(trange+1, 1, ID_front_NI, i), 0, [0.64 0.08 0.18]);
    ploterrorbar(cpawl_all_NI(trange+1, 1, ID_front_NI, i), cpawl_all_NI(trange+1, 3, ID_front_NI, i), [0.00,0.45,0.74]);
    ploterrorbar(cpawr_all_NI(trange+1, 1, ID_front_NI, i), cpawr_all_NI(trange+1, 3, ID_front_NI, i), [0.47,0.67,0.19]);
    ploterrorbar(ctop_all_NI(trange+1, 1, ID_front_NI, i), ctop_all_NI(trange+1, 3, ID_front_NI, i), [0.93,0.69,0.13]);
    ploterrorbar(cbottom_all_NI(trange+1, 1, ID_front_NI, i), cbottom_all_NI(trange+1, 3, ID_front_NI, i), [0.93,0.69,0.13]);
end

% inhibition data
if ~isempty(ID_front_I)
    plotline(0, ceyel_all_I(trange+1, 1, ID_front_I, i), 0, ceyel_all_I(trange+1, 3, ID_front_I, i), [0.8 0.8 0.8]);
    plotline(cnose_all_I(trange+1, 1, ID_front_I, i), ceyel_all_I(trange+1, 1, ID_front_I, i), cnose_all_I(trange+1, 3, ID_front_I, i), ceyel_all_I(trange+1, 3, ID_front_I, i), [0.8 0.8 0.8]);
    plotline(cnose_all_I(trange+1, 1, ID_front_I, i), 0, cnose_all_I(trange+1, 3, ID_front_I, i), 0, [0.8 0.8 0.8]);
    
    plotline(cnose_all_I(trange+1, 1, ID_front_I, i), cpawl_all_I(trange+1, 1, ID_front_I, i), cnose_all_I(trange+1, 3, ID_front_I, i), cpawl_all_I(trange+1, 3, ID_front_I, i), [0.8 0.8 0.8]);
    plotline(cnose_all_I(trange+1, 1, ID_front_I, i), cpawr_all_I(trange+1, 1, ID_front_I, i), cnose_all_I(trange+1, 3, ID_front_I, i), cpawr_all_I(trange+1, 3, ID_front_I, i), [0.8 0.8 0.8]);
    plotline(cpawl_all_I(trange+1, 1, ID_front_I, i), cpawr_all_I(trange+1, 1, ID_front_I, i), cpawl_all_I(trange+1, 3, ID_front_I, i), cpawr_all_I(trange+1, 3, ID_front_I, i), [0.8 0.8 0.8]);
    
    plotline(ctop_all_I(trange+1, 1, ID_front_I, i), cbottom_all_I(trange+1, 1, ID_front_I, i),...
        ctop_all_I(trange+1, 3, ID_front_I, i), cbottom_all_I(trange+1, 3, ID_front_I, i), [0.92,0.92,0.45]);
    
    ploterrorbar(ceyel_all_I(trange+1, 1, ID_front_I, i), 0, [0.8 0.8 0.8]);
    ploterrorbar(cnose_all_I(trange+1, 1, ID_front_I, i), 0, [1 0 0]);
    ploterrorbar(cpawl_all_I(trange+1, 1, ID_front_I, i), cpawl_all_I(trange+1, 3, ID_front_I, i), [0.30,0.75,0.93]);
    ploterrorbar(cpawr_all_I(trange+1, 1, ID_front_I, i), cpawr_all_I(trange+1, 3, ID_front_I, i), [0 1 0]);
    ploterrorbar(ctop_all_I(trange+1, 1, ID_front_I, i), ctop_all_I(trange+1, 3, ID_front_I, i), [0.92,0.92,0.45]);
    ploterrorbar(cbottom_all_I(trange+1, 1, ID_front_I, i), cbottom_all_I(trange+1, 3, ID_front_I, i), [0.92,0.92,0.45]);
end

axis equal;
xlabel('X (mm)');
ylabel('Z (mm)');
set(gca, 'TickDir', 'out', 'FontSize', 12, 'ButtonDownFcn', @extract_figure);

subplot(1, 2, 2);
hold on;
if ~isempty(ID_top_NI)
    plotline(0, ceyel_all_NI(trange+1, 1, ID_top_NI, i), 0, ceyel_all_NI(trange+1, 2, ID_top_NI, i), [0 0 0]);
    plotline(cnose_all_NI(trange+1, 1, ID_top_NI, i), ceyel_all_NI(trange+1, 1, ID_top_NI, i), cnose_all_NI(trange+1, 2, ID_top_NI, i), ceyel_all_NI(trange+1, 2, ID_top_NI, i), [0 0 0]);
    plotline(cnose_all_NI(trange+1, 1, ID_top_NI, i), 0, cnose_all_NI(trange+1, 2, ID_top_NI, i), 0, [0 0 0]);
    
    plotline(cnose_all_NI(trange+1, 1, ID_top_NI, i), cpawl_all_NI(trange+1, 1, ID_top_NI, i), cnose_all_NI(trange+1, 2, ID_top_NI, i), cpawl_all_NI(trange+1, 2, ID_top_NI, i), [0 0 0]);
    plotline(cnose_all_NI(trange+1, 1, ID_top_NI, i), cpawr_all_NI(trange+1, 1, ID_top_NI, i), cnose_all_NI(trange+1, 2, ID_top_NI, i), cpawr_all_NI(trange+1, 2, ID_top_NI, i), [0 0 0]);
    plotline(cpawl_all_NI(trange+1, 1, ID_top_NI, i), cpawr_all_NI(trange+1, 1, ID_top_NI, i), cpawl_all_NI(trange+1, 2, ID_top_NI, i), cpawr_all_NI(trange+1, 2, ID_top_NI, i), [0 0 0]);
    
    plotline(ctop_all_NI(trange+1, 1, ID_top_NI, i), cbottom_all_NI(trange+1, 1, ID_top_NI, i),...
        ctop_all_NI(trange+1, 2, ID_top_NI, i), cbottom_all_NI(trange+1, 2, ID_top_NI, i), [0.93,0.69,0.13]);

    ploterrorbar(ceyel_all_NI(trange+1, 1, ID_top_NI, i), 0, [0 0 0]);
    ploterrorbar(cnose_all_NI(trange+1, 1, ID_top_NI, i), cnose_all_NI(trange+1, 2, ID_top_NI, i), [0.64 0.08 0.18]);
    ploterrorbar(cpawl_all_NI(trange+1, 1, ID_top_NI, i), cpawl_all_NI(trange+1, 2, ID_top_NI, i), [0.00,0.45,0.74]);
    ploterrorbar(cpawr_all_NI(trange+1, 1, ID_top_NI, i), cpawr_all_NI(trange+1, 2, ID_top_NI, i), [0.47,0.67,0.19]);
    ploterrorbar(ctop_all_NI(trange+1, 1, ID_top_NI, i), ctop_all_NI(trange+1, 2, ID_top_NI, i), [0.93,0.69,0.13]);
    ploterrorbar(cbottom_all_NI(trange+1, 1, ID_top_NI, i), cbottom_all_NI(trange+1, 2, ID_top_NI, i), [0.93,0.69,0.13]);
end

% inhibition data
if ~isempty(ID_top_I)
    plotline(0, ceyel_all_I(trange+1, 1, ID_top_I, i), 0, ceyel_all_I(trange+1, 2, ID_top_I, i), [0.8 0.8 0.8]);
    plotline(cnose_all_I(trange+1, 1, ID_top_I, i), ceyel_all_I(trange+1, 1, ID_top_I, i), cnose_all_I(trange+1, 2, ID_top_I, i), ceyel_all_I(trange+1, 2, ID_top_I, i), [0.8 0.8 0.8]);
    plotline(cnose_all_I(trange+1, 1, ID_top_I, i), 0, cnose_all_I(trange+1, 2, ID_top_I, i), 0, [0.8 0.8 0.8]);
    
    plotline(cnose_all_I(trange+1, 1, ID_top_I, i), cpawl_all_I(trange+1, 1, ID_top_I, i), cnose_all_I(trange+1, 2, ID_top_I, i), cpawl_all_I(trange+1, 2, ID_top_I, i), [0.8 0.8 0.8]);
    plotline(cnose_all_I(trange+1, 1, ID_top_I, i), cpawr_all_I(trange+1, 1, ID_top_I, i), cnose_all_I(trange+1, 2, ID_top_I, i), cpawr_all_I(trange+1, 2, ID_top_I, i), [0.8 0.8 0.8]);
    plotline(cpawl_all_I(trange+1, 1, ID_top_I, i), cpawr_all_I(trange+1, 1, ID_top_I, i), cpawl_all_I(trange+1, 2, ID_top_I, i), cpawr_all_I(trange+1, 2, ID_top_I, i), [0.8 0.8 0.8]);
    
    plotline(ctop_all_I(trange+1, 1, ID_top_I, i), cbottom_all_I(trange+1, 1, ID_top_I, i),...
        ctop_all_I(trange+1, 2, ID_top_I, i), cbottom_all_I(trange+1, 2, ID_top_I, i), [0.92,0.92,0.45]);
    
    ploterrorbar(ceyel_all_I(trange+1, 1, ID_top_I, i), 0, [0.8 0.8 0.8]);
    ploterrorbar(cnose_all_I(trange+1, 1, ID_top_I, i), cnose_all_I(trange+1, 2, ID_top_I, i), [1 0 0]);
    ploterrorbar(cpawl_all_I(trange+1, 1, ID_top_I, i), cpawl_all_I(trange+1, 2, ID_top_I, i), [0.30,0.75,0.93]);
    ploterrorbar(cpawr_all_I(trange+1, 1, ID_top_I, i), cpawr_all_I(trange+1, 2, ID_top_I, i), [0 1 0]);
    ploterrorbar(ctop_all_I(trange+1, 1, ID_top_I, i), ctop_all_I(trange+1, 2, ID_top_I, i), [0.92,0.92,0.45]);
    ploterrorbar(cbottom_all_I(trange+1, 1, ID_top_I, i), cbottom_all_I(trange+1, 2, ID_top_I, i), [0.92,0.92,0.45]);
end

axis equal;
xlabel('X (mm)');
ylabel('Y (mm)');
set(gca, 'TickDir', 'out', 'FontSize', 12, 'ButtonDownFcn', @extract_figure);