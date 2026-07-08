function articleFormatMultiplePlots(fontName, subH)
    % ARTICLEFORMATMULTIPLEPLOTS Professional Academic Version
    %   - Each subplot: 3.3" width x 2.0" height.
    %   - PROTECTS: Externally defined Y-limits (ylim).
    %   - TIGHTENS: Only the X-axis (Capacity/SOC).
    %   - FIXES: Hidden yyaxis right labels with smart padding.

    if ~exist("fontName", "var"); fontName = 'Times New Roman'; end
    if ~exist("subH", "var"); subH = 3.5; end

    % --- 1. Capture figure and axes ---
    fig = gcf;
    allAxes = findobj(fig, 'Type', 'axes', '-not', 'Tag', 'legend', '-not', 'Tag', 'Colorbar');
    if isempty(allAxes), return; end

    % --- 2. Dynamic Figure Sizing ---
    subW = 8; 

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

    % --- 3. TiledLayout Spacing ---
    tLayout = findobj(fig, 'Type', 'tiledlayout');
    if ~isempty(tLayout)
        tLayout.TileSpacing = 'none';
        tLayout.Padding = 'none';
    end

    % --- 4. Loop through each subplot ---
    for i = 1:length(allAxes)
        ax = allAxes(i);
        
        % --- STEP A: Protect your Y-limits before tightening ---
        if length(ax.YAxis) > 1
            % Save Left and Right limits for dual-axes
            limL = ax.YAxis(1).Limits;
            limR = ax.YAxis(2).Limits;
            xlim(ax, 'tight'); % Only tighten X
            ax.YAxis(1).Limits = limL; % Restore
            ax.YAxis(2).Limits = limR; % Restore
        else
            limY = ax.YLim;
            xlim(ax, 'tight'); % Only tighten X
            ax.YLim = limY;    % Restore
        end
        
        % --- STEP B: YYAXIS Support & Formatting ---
        for ySide = 1:length(ax.YAxis)
            yObj = ax.YAxis(ySide);
            % yObj.Color = [0 0 0]; 
            yObj.Label.Interpreter = 'latex';
            yObj.FontSize = 8;
            yObj.Label.FontName = fontName;
        end

        % --- STEP C: X-Axis Styling ---
        ax.XAxis.Color = [0 0 0];
        ax.XAxis.FontSize = 8;
        ax.XLabel.FontName = fontName;
        ax.XLabel.Interpreter = 'latex';

        % --- STEP D: General Axes Styling ---
        title(ax, ''); 
        ax.Box = 'on';           
        ax.LineWidth = 0.7;      
        ax.FontName = fontName;
        ax.XMinorTick = 'on';
        ax.YMinorTick = 'on';
        ax.TickLabelInterpreter = 'latex';
        
        % --- STEP E: Data Line Styling ---
        lines = findobj(ax, 'Type', 'line');
        for k = 1:length(lines)
            if lines(k).LineWidth ==  0.75
                lines(k).LineWidth = 2;
            end
        end
        % set(lines(lines.LineWidth == 1), 'LineWidth', 2);
        % set(findobj(ax, 'Type', 'line'), 'LineWidth', 2); 
        
        % --- STEP F: Grid ---
        grid(ax, 'on');
        ax.GridLineStyle = ':';
        ax.GridAlpha = 0.15;

        % % --- STEP G: SMART CROP (Fix for hidden yyaxis labels) ---
        % drawnow; % Required for accurate TightInset
        % ti = get(ax, 'TightInset');
        % 
        % if length(ax.YAxis) > 1
        %     % Added padding specifically for the right-side label
        %     extraRightPadding = 2; 
        %     set(ax, 'LooseInset', [ti(1), ti(2), ti(3) + extraRightPadding, ti(4)]);
        % else
        %     set(ax, 'LooseInset', ti);
        % end
    end

    drawnow;
end