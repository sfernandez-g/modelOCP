%PPTFORMATMULTIPLEPLOTS Format the current figure for PowerPoint export.
%
%   pptFormatMultiplePlots()
%   pptFormatMultiplePlots(fontName)
%   pptFormatMultiplePlots(fontName, subH)
%
%   Applies a consistent PowerPoint-oriented format to all axes in the
%   current figure, including figures containing tiled layouts.
%
%   Inputs
%   -------
%   fontName : char/string, optional
%       Font used for axes and labels.
%       Default: 'Times New Roman'.
%
%   subH : scalar, optional
%       Height assigned to each subplot in centimeters.
%       Default: 7 cm.
%
%   Notes
%   -----
%   The subplot width is fixed to 16 cm. The total figure size is inferred
%   from the tiled-layout dimensions or, for standard axes, from their
%   positions.
%
%   The function also:
%       - sets figure and paper dimensions in centimeters;
%       - removes tiled-layout spacing and padding;
%       - preserves existing Y limits while tightening X limits;
%       - applies LaTeX tick and label interpreters;
%       - enables minor ticks and a light dotted grid;
%       - increases 1-pt line widths to 2 pt.

function pptFormatMultiplePlots(fontName, subH)

    if ~exist("fontName", "var") || isempty(fontName)
        fontName = 'Times New Roman';
    end

    if ~exist("subH", "var") || isempty(subH)
        subH = 7;
    end

    fig = gcf;

    allAxes = findobj(fig, ...
        'Type', 'axes', ...
        '-not', 'Tag', 'legend', ...
        '-not', 'Tag', 'Colorbar');

    if isempty(allAxes)
        return;
    end

    subW = 16;

    % Determine figure grid dimensions.
    tLayout = findobj(fig, 'Type', 'tiledlayout');

    if ~isempty(tLayout)
        tLayout = tLayout(1);

        tLayout.TileSpacing = 'none';
        tLayout.Padding = 'none';

        gridSize = tLayout.GridSize;
        numRows = gridSize(1);
        numCols = gridSize(2);
    else
        allPos = get(allAxes, 'Position');

        if iscell(allPos)
            posMatrix = vertcat(allPos{:});
        else
            posMatrix = allPos;
        end

        numCols = numel(unique(round(posMatrix(:, 1), 2)));
        numRows = numel(unique(round(posMatrix(:, 2), 2)));
    end

    totalW = subW * numCols;
    totalH = subH * numRows;

    % Figure and export properties.
    fig.Units = 'centimeters';
    fig.Position = [10 10 totalW totalH];
    fig.Color = [1 1 1];

    fig.PaperUnits = 'centimeters';
    fig.PaperPositionMode = 'auto';
    fig.PaperSize = [totalW totalH];

    % Format axes.
    for i = 1:numel(allAxes)

        ax = allAxes(i);

        hasBar = ~isempty(findobj(ax, 'Type', 'bar'));
        isCat = iscategorical(ax.XAxis);

        % Tighten X limits without changing Y limits.
        if ~hasBar && ~isCat

            if numel(ax.YAxis) > 1

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

        % Axis labels and fonts.
        for ySide = 1:numel(ax.YAxis)
            yObj = ax.YAxis(ySide);
            yObj.Label.Interpreter = 'latex';
            yObj.Label.FontName = fontName;
            yObj.FontSize = 16;
        end

        ax.XAxis.Color = [0 0 0];
        ax.XAxis.FontSize = 16;
        ax.XLabel.FontName = fontName;
        ax.XLabel.Interpreter = 'latex';

        ax.FontName = fontName;
        ax.TickLabelInterpreter = 'latex';

        % Axes appearance.
        title(ax, '');
        ax.Box = 'on';
        ax.LineWidth = 0.7;

        ax.XMinorTick = 'on';
        ax.YMinorTick = 'on';

        % Line widths.
        lines = findobj(ax, 'Type', 'line');

        for k = 1:numel(lines)
            if lines(k).LineWidth == 1
                lines(k).LineWidth = 2;
            end
        end

        % Grid.
        grid(ax, 'on');
        ax.GridLineStyle = ':';
        ax.GridAlpha = 0.15;

    end

    drawnow;

end