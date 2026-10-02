function pptFormatMultiplePlots(fontName, subH)
    if ~exist("fontName", "var"); fontName = 'Times New Roman'; end
    if ~exist("subH", "var"); subH = 7; end
    
    fig = gcf;
    allAxes = findobj(fig, 'Type', 'axes', '-not', 'Tag', 'legend', '-not', 'Tag', 'Colorbar');
    if isempty(allAxes), return; end
    subW = 16; 
    allPos = get(allAxes, 'Position');
    if iscell(allPos), posMatrix = vertcat(allPos{:}); else, posMatrix = allPos; end
    
    numCols = length(unique(round(posMatrix(:,1), 2))); 
    numRows = length(unique(round(posMatrix(:,2), 2))); 
    totalW = subW * numCols;
    totalH = subH * numRows;
    fig.Units = 'centimeters';
    fig.Position = [10 10 totalW totalH]; 
    fig.Color = [1 1 1];        
    fig.PaperUnits = 'centimeters';
    fig.PaperPositionMode = 'auto'; 
    fig.PaperSize = [totalW totalH]; 
    tLayout = findobj(fig, 'Type', 'tiledlayout');
    if ~isempty(tLayout)
        tLayout.TileSpacing = 'none';
        tLayout.Padding = 'none';
    end
    for i = 1:length(allAxes)
        ax = allAxes(i);
        
        hasBar = ~isempty(findobj(ax, 'Type', 'bar'));
        isCat  = iscategorical(ax.XAxis);
        
        if ~hasBar && ~isCat
            if length(ax.YAxis) > 1
                limL = ax.YAxis(1).Limits;
                limR = ax.YAxis(2).Limits;
                xlim(ax, 'tight');
                ax.YAxis(1).Limits = limL;
                ax.YAxis(2).Limits = limR;
            else
                limY = ax.YLim;
                xlim(ax, 'tight');
                ax.YLim = limY;
            end
        end
        
        for ySide = 1:length(ax.YAxis)
            yObj = ax.YAxis(ySide);
            yObj.Label.Interpreter = 'latex';
            yObj.FontSize = 16;
            yObj.Label.FontName = fontName;
        end
        ax.XAxis.Color = [0 0 0];
        ax.XAxis.FontSize = 16;
        ax.XLabel.FontName = fontName;
        ax.XLabel.Interpreter = 'latex';
        title(ax, ''); 
        ax.Box = 'on';           
        ax.LineWidth = 0.7;      
        ax.FontName = fontName;
        ax.XMinorTick = 'on';
        ax.YMinorTick = 'on';
        ax.TickLabelInterpreter = 'latex';
        
        if ~isempty(ax.XTickLabel)
            lbls = cellstr(ax.XTickLabel);
            for j = 1:length(lbls)
                if (contains(lbls{j}, '_') || contains(lbls{j}, '^')) && ~contains(lbls{j}, '$')
                    lbls{j} = ['$' lbls{j} '$'];
                end
            end
            ax.XTickLabel = lbls;
        end
        lines = findobj(ax, 'Type', 'line');
        for k = 1:length(lines)
            if lines(k).LineWidth == 1
                lines(k).LineWidth = 2;
            end
        end
        
        grid(ax, 'on');
        ax.GridLineStyle = ':';
        ax.GridAlpha = 0.15;
    end
    drawnow;
end