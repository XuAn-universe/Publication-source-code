function draw_individual_locations_at_bite(ceyel_all_NI, ceyer_all_NI, cpawl_all_NI, cpawr_all_NI, ctop_all_NI, cbottom_all_NI,...
    ceyel_all_I, ceyer_all_I, cpawl_all_I, cpawr_all_I, ctop_all_I, cbottom_all_I, ID_front_NI, ID_top_NI, ID_front_I, ID_top_I, trange, i, biteinfo_NI, biteinfo_I)
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
    for j = ID_front_NI
        plot([ceyer_all_NI(trange+1, 1, j, i) ceyel_all_NI(trange+1, 1, j, i)], [ceyer_all_NI(trange+1, 3, j, i) ceyel_all_NI(trange+1, 3, j, i)], '-k',...
            'Marker', '.', 'MarkerSize', 12, 'ButtonDownFcn', @displayfilename, 'UserData', biteinfo_NI{j});
        plot([0 ceyel_all_NI(trange+1, 1, j, i)], [0 ceyel_all_NI(trange+1, 3, j, i)], '-k', 'ButtonDownFcn', @displayfilename, 'UserData', biteinfo_NI{j});
        plot([0 ceyer_all_NI(trange+1, 1, j, i)], [0 ceyer_all_NI(trange+1, 3, j, i)], '-k', 'ButtonDownFcn', @displayfilename, 'UserData', biteinfo_NI{j});
        
        if ~any(isnan([cpawl_all_NI(trange+1, 1, j, i) cpawl_all_NI(trange+1, 3, j, i) cpawr_all_NI(trange+1, 1, j, i) cpawr_all_NI(trange+1, 3, j, i)]))
            plot([0 cpawl_all_NI(trange+1, 1, j, i)], [0 cpawl_all_NI(trange+1, 3, j, i)], '-k', 'ButtonDownFcn', @displayfilename, 'UserData', biteinfo_NI{j});
            plot([0 cpawr_all_NI(trange+1, 1, j, i)], [0 cpawr_all_NI(trange+1, 3, j, i)], '-k', 'ButtonDownFcn', @displayfilename, 'UserData', biteinfo_NI{j});
            plot([cpawl_all_NI(trange+1, 1, j, i) cpawr_all_NI(trange+1, 1, j, i)], [cpawl_all_NI(trange+1, 3, j, i) cpawr_all_NI(trange+1, 3, j, i)], '-k', 'ButtonDownFcn', @displayfilename, 'UserData', biteinfo_NI{j});
            
            plot(cpawl_all_NI(trange+1, 1, j, i), cpawl_all_NI(trange+1, 3, j, i), '.', 'Color', [0.00,0.45,0.74], 'MarkerSize', 12);
            plot(cpawr_all_NI(trange+1, 1, j, i), cpawr_all_NI(trange+1, 3, j, i), '.', 'Color', [0.47,0.67,0.19], 'MarkerSize', 12);
        end
        
        if ~any(isnan([ctop_all_NI(trange+1, 1, j, i) ctop_all_NI(trange+1, 3, j, i) cbottom_all_NI(trange+1, 1, j, i) cbottom_all_NI(trange+1, 3, j, i)]))
            plot([ctop_all_NI(trange+1, 1, j, i) cbottom_all_NI(trange+1, 1, j, i)], [ctop_all_NI(trange+1, 3, j, i) cbottom_all_NI(trange+1, 3, j, i)], '-',...
                'Marker', '.', 'Color', [0.93,0.69,0.13], 'MarkerSize', 12, 'ButtonDownFcn', @displayfilename, 'UserData', biteinfo_NI{j});
        end
    end
end
if ~isempty(ID_front_I)
    for j = ID_front_I
        plot([ceyer_all_I(trange+1, 1, j, i) ceyel_all_I(trange+1, 1, j, i)], [ceyer_all_I(trange+1, 3, j, i) ceyel_all_I(trange+1, 3, j, i)], '-',...
            'Marker', '.', 'Color', [0.8 0.8 0.8], 'MarkerSize', 12, 'ButtonDownFcn', @displayfilename, 'UserData', biteinfo_I{j});
        plot([0 ceyel_all_I(trange+1, 1, j, i)], [0 ceyel_all_I(trange+1, 3, j, i)], '-', 'Color', [0.8 0.8 0.8], 'ButtonDownFcn', @displayfilename, 'UserData', biteinfo_I{j});
        plot([0 ceyer_all_I(trange+1, 1, j, i)], [0 ceyer_all_I(trange+1, 3, j, i)], '-', 'Color', [0.8 0.8 0.8], 'ButtonDownFcn', @displayfilename, 'UserData', biteinfo_I{j});
        
        if ~any(isnan([cpawl_all_I(trange+1, 1, j, i) cpawl_all_I(trange+1, 3, j, i) cpawr_all_I(trange+1, 1, j, i) cpawr_all_I(trange+1, 3, j, i)]))
            plot([0 cpawl_all_I(trange+1, 1, j, i)], [0 cpawl_all_I(trange+1, 3, j, i)], '-', 'Color', [0.8 0.8 0.8], 'ButtonDownFcn', @displayfilename, 'UserData', biteinfo_I{j});
            plot([0 cpawr_all_I(trange+1, 1, j, i)], [0 cpawr_all_I(trange+1, 3, j, i)], '-', 'Color', [0.8 0.8 0.8], 'ButtonDownFcn', @displayfilename, 'UserData', biteinfo_I{j});
            plot([cpawl_all_I(trange+1, 1, j, i) cpawr_all_I(trange+1, 1, j, i)], [cpawl_all_I(trange+1, 3, j, i) cpawr_all_I(trange+1, 3, j, i)], '-',...
                'Color', [0.8 0.8 0.8], 'ButtonDownFcn', @displayfilename, 'UserData', biteinfo_I{j});
            
            plot(cpawl_all_I(trange+1, 1, j, i), cpawl_all_I(trange+1, 3, j, i), '.', 'Color', [0.30,0.75,0.93], 'MarkerSize', 12);
            plot(cpawr_all_I(trange+1, 1, j, i), cpawr_all_I(trange+1, 3, j, i), '.g', 'MarkerSize', 12);
        end
        
        if ~any(isnan([ctop_all_I(trange+1, 1, j, i) ctop_all_I(trange+1, 3, j, i) cbottom_all_I(trange+1, 1, j, i) cbottom_all_I(trange+1, 3, j, i)]))
            plot([ctop_all_I(trange+1, 1, j, i) cbottom_all_I(trange+1, 1, j, i)], [ctop_all_I(trange+1, 3, j, i) cbottom_all_I(trange+1, 3, j, i)], '-',...
                'Marker', '.', 'Color', [0.92,0.92,0.45], 'MarkerSize', 12, 'ButtonDownFcn', @displayfilename, 'UserData', biteinfo_I{j});
        end
    end
end
axis equal;
xlabel('X (mm)');
ylabel('Z (mm)');
set(gca, 'TickDir', 'out', 'FontSize', 12, 'ButtonDownFcn', @extract_figure);

subplot(1, 2, 2);
hold on;
if ~isempty(ID_top_NI)
    for j = ID_top_NI
        plot([ceyer_all_NI(trange+1, 1, j, i) ceyel_all_NI(trange+1, 1, j, i)], [ceyer_all_NI(trange+1, 2, j, i) ceyel_all_NI(trange+1, 2, j, i)], '-k',...
            'Marker', '.', 'MarkerSize', 12, 'ButtonDownFcn', @displayfilename, 'UserData', biteinfo_NI{j});
        plot([0 ceyel_all_NI(trange+1, 1, j, i)], [0 ceyel_all_NI(trange+1, 2, j, i)], '-k', 'ButtonDownFcn', @displayfilename, 'UserData', biteinfo_NI{j});
        plot([0 ceyer_all_NI(trange+1, 1, j, i)], [0 ceyer_all_NI(trange+1, 2, j, i)], '-k', 'ButtonDownFcn', @displayfilename, 'UserData', biteinfo_NI{j});
        
        if ~any(isnan([cpawl_all_NI(trange+1, 1, j, i) cpawl_all_NI(trange+1, 2, j, i) cpawr_all_NI(trange+1, 1, j, i) cpawr_all_NI(trange+1, 2, j, i)]))
            plot([0 cpawl_all_NI(trange+1, 1, j, i)], [0 cpawl_all_NI(trange+1, 2, j, i)], '-k', 'ButtonDownFcn', @displayfilename, 'UserData', biteinfo_NI{j});
            plot([0 cpawr_all_NI(trange+1, 1, j, i)], [0 cpawr_all_NI(trange+1, 2, j, i)], '-k', 'ButtonDownFcn', @displayfilename, 'UserData', biteinfo_NI{j});
            plot([cpawl_all_NI(trange+1, 1, j, i) cpawr_all_NI(trange+1, 1, j, i)], [cpawl_all_NI(trange+1, 2, j, i) cpawr_all_NI(trange+1, 2, j, i)], '-k', 'ButtonDownFcn', @displayfilename, 'UserData', biteinfo_NI{j});
            
            plot(cpawl_all_NI(trange+1, 1, j, i), cpawl_all_NI(trange+1, 2, j, i), '.', 'Color', [0.00,0.45,0.74], 'MarkerSize', 12);
            plot(cpawr_all_NI(trange+1, 1, j, i), cpawr_all_NI(trange+1, 2, j, i), '.', 'Color', [0.47,0.67,0.19], 'MarkerSize', 12);
        end
        
        if ~any(isnan([ctop_all_NI(trange+1, 1, j, i) ctop_all_NI(trange+1, 2, j, i) cbottom_all_NI(trange+1, 1, j, i) cbottom_all_NI(trange+1, 2, j, i)]))
            plot([ctop_all_NI(trange+1, 1, j, i) cbottom_all_NI(trange+1, 1, j, i)], [ctop_all_NI(trange+1, 2, j, i) cbottom_all_NI(trange+1, 2, j, i)], '-',...
                'Marker', '.', 'Color', [0.93,0.69,0.13], 'MarkerSize', 12, 'ButtonDownFcn', @displayfilename, 'UserData', biteinfo_NI{j});
        end
    end
end
if ~isempty(ID_top_I)
    for j = ID_top_I
        plot([ceyer_all_I(trange+1, 1, j, i) ceyel_all_I(trange+1, 1, j, i)], [ceyer_all_I(trange+1, 2, j, i) ceyel_all_I(trange+1, 2, j, i)], '-',...
            'Marker', '.', 'Color', [0.8 0.8 0.8], 'MarkerSize', 12, 'ButtonDownFcn', @displayfilename, 'UserData', biteinfo_I{j});
        plot([0 ceyel_all_I(trange+1, 1, j, i)], [0 ceyel_all_I(trange+1, 2, j, i)], '-', 'Color', [0.8 0.8 0.8], 'ButtonDownFcn', @displayfilename, 'UserData', biteinfo_I{j});
        plot([0 ceyer_all_I(trange+1, 1, j, i)], [0 ceyer_all_I(trange+1, 2, j, i)], '-', 'Color', [0.8 0.8 0.8], 'ButtonDownFcn', @displayfilename, 'UserData', biteinfo_I{j});
        
        if ~any(isnan([cpawl_all_I(trange+1, 1, j, i) cpawl_all_I(trange+1, 2, j, i) cpawr_all_I(trange+1, 1, j, i) cpawr_all_I(trange+1, 2, j, i)]))
            plot([0 cpawl_all_I(trange+1, 1, j, i)], [0 cpawl_all_I(trange+1, 2, j, i)], '-', 'Color', [0.8 0.8 0.8], 'ButtonDownFcn', @displayfilename, 'UserData', biteinfo_I{j});
            plot([0 cpawr_all_I(trange+1, 1, j, i)], [0 cpawr_all_I(trange+1, 2, j, i)], '-', 'Color', [0.8 0.8 0.8], 'ButtonDownFcn', @displayfilename, 'UserData', biteinfo_I{j});
            plot([cpawl_all_I(trange+1, 1, j, i) cpawr_all_I(trange+1, 1, j, i)], [cpawl_all_I(trange+1, 2, j, i) cpawr_all_I(trange+1, 2, j, i)], '-',...
                'Color', [0.8 0.8 0.8], 'ButtonDownFcn', @displayfilename, 'UserData', biteinfo_I{j});
            
            plot(cpawl_all_I(trange+1, 1, j, i), cpawl_all_I(trange+1, 2, j, i), '.', 'Color', [0.30,0.75,0.93], 'MarkerSize', 12);
            plot(cpawr_all_I(trange+1, 1, j, i), cpawr_all_I(trange+1, 2, j, i), '.g', 'MarkerSize', 12);
        end
        
        if ~any(isnan([ctop_all_I(trange+1, 1, j, i) ctop_all_I(trange+1, 2, j, i) cbottom_all_I(trange+1, 1, j, i) cbottom_all_I(trange+1, 2, j, i)]))
            plot([ctop_all_I(trange+1, 1, j, i) cbottom_all_I(trange+1, 1, j, i)], [ctop_all_I(trange+1, 2, j, i) cbottom_all_I(trange+1, 2, j, i)], '-',...
                'Marker', '.', 'Color', [0.92,0.92,0.45], 'MarkerSize', 12, 'ButtonDownFcn', @displayfilename, 'UserData', biteinfo_I{j});
        end
    end
end
axis equal;
xlabel('X (mm)');
ylabel('Y (mm)');
set(gca, 'TickDir', 'out', 'FontSize', 12, 'ButtonDownFcn', @extract_figure);

    function displayfilename(src, eventdata)
        htext = text(eventdata.IntersectionPoint(1), eventdata.IntersectionPoint(2), src.UserData, 'FontSize', 12, 'Color', 'r');
        pause(2);
        try
            delete(htext)
            clear htext;
        end
    end
end