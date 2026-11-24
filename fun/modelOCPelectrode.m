function [model, ocvrel, supplementaryResults, ocpFigures] = modelOCPelectrode(settings)
%MODELOCPELECTRODE Process half-cell OCP data and fit an electrode MSMR model.
%   [model, ocvrel, supplementaryResults, ocpFigures] = modelOCPelectrode(settings)
%
%   settings is a scalar structure with electrode-specific fields:
%     inputFile              Arbin XLSX path.
%     capacityFile           MAT path containing scalar Q_r [mAh].
%     cellName               Half-cell identifier passed to processing.
%     electrodeName          Electrode label used in figure window names.
%     voltageRange           [min max] potential versus Li [V].
%     activeMaterialFraction Initial active material mass fraction (0, 1].
%     lithiationAnchor       "delithiated" (theta_min=0, charge branch) or
%                            "lithiated" (theta_max=1, discharge branch).
%     initialMSMR            Struct with U0 [V], X [-], w, or [] to guess.
%     galleryCount           Number of galleries when initialMSMR is [].
%     boundScales            [U0 X w] multipliers for lower bounds;
%                            upper bounds use their reciprocals.
%     differentialWeight     Weight of the differential residual.
%     hiddenGalleryCorrection Logical enabling the original NMC correction.
%   Optional fields:
%     temperatureDegC        Temperature [degC], default 25.
%     relativeInsets, fitInsets, absoluteInsets
%                            N-by-6 rows [left bottom width height xmin xmax].
%                            Positions normalized; lithiation limits [-].
%                            Default zeros(0,6), disabling insets.
%     showRelativeLegend     Figure 2 legend, default true.
%     showComparisonLegend   Figure 7 legend, default true.
%     outputFile             MAT output path; default "" disables saving.
%     optimizer              Settings forwarded to fitMSMR.
%
%   Outputs preserve the original model/LUT schema. ocvrel contains processed
%   data and fitted curves; supplementaryResults snapshots results before
%   corrections; ocpFigures is a 7-by-1 array of figure handles numbered
%   consecutively. Figure 6 contains three MSMR composition/gallery subplots.
%   See example_modelOCPelectrode_graphite and example_modelOCPelectrode_NMC
%   for complete configurations.
%   This function creates figures but does not clear the caller workspace.

    arguments
        settings (1,1) struct
    end
    required = {'inputFile', 'capacityFile', 'cellName', 'electrodeName', ...
        'voltageRange', 'activeMaterialFraction', 'lithiationAnchor', ...
        'initialMSMR', 'boundScales', 'differentialWeight', ...
        'hiddenGalleryCorrection'};
    for i = 1:numel(required)
        if ~isfield(settings, required{i})
            error('modelOCPelectrode:MissingSetting', ...
                'Missing settings.%s.', required{i});
        end
    end
    defaults = struct('temperatureDegC', 25, 'relativeInsets', zeros(0,6), ...
        'fitInsets', zeros(0,6), 'absoluteInsets', zeros(0,6), ...
        'showRelativeLegend', true, 'showComparisonLegend', true, ...
        'outputFile', "");
    names = fieldnames(defaults);
    for i = 1:numel(names)
        if ~isfield(settings, names{i})
            settings.(names{i}) = defaults.(names{i});
        end
    end
    settings.lithiationAnchor = string(settings.lithiationAnchor);
    mustBeMember(settings.lithiationAnchor, ["delithiated", "lithiated"]);
    validateattributes(settings.lithiationAnchor, {'string'}, {'scalar'});
    validateattributes(settings.voltageRange, {'numeric'}, ...
        {'vector', 'numel', 2, 'real', 'finite', 'increasing'});
    validateattributes(settings.activeMaterialFraction, {'numeric'}, ...
        {'scalar', 'real', 'finite', 'positive', '<=', 1});
    validateattributes(settings.boundScales, {'numeric'}, ...
        {'vector', 'numel', 3, 'real', 'finite', 'positive', '<=', 1});
    validateattributes(settings.temperatureDegC, {'numeric'}, ...
        {'scalar', 'real', 'finite', '>', -273.15});
    validateattributes(settings.differentialWeight, {'numeric'}, ...
        {'scalar', 'real', 'finite', 'nonnegative'});
    for name = ["hiddenGalleryCorrection", "showRelativeLegend", "showComparisonLegend"]
        validateattributes(settings.(name), {'logical', 'numeric'}, ...
            {'scalar', 'real', 'finite', 'binary'});
    end
    if isempty(settings.initialMSMR)
        if ~isfield(settings, 'galleryCount')
            error('modelOCPelectrode:MissingGalleryCount', ...
                'settings.galleryCount is required when initialMSMR is empty.');
        end
        validateattributes(settings.galleryCount, {'numeric'}, ...
            {'scalar', 'integer', 'positive', 'finite'});
    end
    for name = ["relativeInsets", "fitInsets", "absoluteInsets"]
        validateattributes(settings.(name), {'numeric'}, ...
            {'2d', 'ncols', 6, 'real', 'finite'});
        insets = settings.(name);
        if any(insets(:,3:4) <= 0, 'all') || any(insets(:,5) >= insets(:,6))
            error('modelOCPelectrode:InvalidInsets', ...
                '%s requires positive sizes and increasing lithiation limits.', name);
        end
    end

    repoRoot = fileparts(fileparts(mfilename('fullpath')));
    ocpFigures = gobjects(7, 1); % Indices match the figure window names.
    temperatureDegC = settings.temperatureDegC;
    gasConstant = 8.3144598; % J/(mol K); same constants as compositionMSMR.
    faradayConstant = 96485.3329; % C/mol.
    thermalFactor = faradayConstant/(gasConstant*(temperatureDegC + 273.15)); % 1/V.

    %% add utility paths
    addpath(fullfile(repoRoot, "fun"));

    %% process OCP data
    file = settings.inputFile;
    %--------------------------------------------------------------------------
    cellName = settings.cellName; 
    %--------------------------------------------------------------------------
    v_range = settings.voltageRange;
    % -------------------------------------------------------------------------
    sol_range   = [0 1];
    % -------------------------------------------------------------------------
    capacity = load(settings.capacityFile, 'Q_r');
    validateattributes(capacity.Q_r, {'numeric'}, {'scalar', 'real', 'finite', 'positive'});
    Q_r = capacity.Q_r*1e-3; % Capacity file stores Q_r in mAh; use Ah below.
    AM_per = settings.activeMaterialFraction;
    Q_r_AM  = Q_r*AM_per;
    % -------------------------------------------------------------------------
    test_data = getDataArbin_XLSX(file, 7); % filename, step where TC increments
    % -------------------------------------------------------------------------
    ocpFigures(1) = figure('Name', sprintf('1. %s: Measured half-cell voltage', settings.electrodeName), 'NumberTitle', 'off');
    plot(test_data.t./3600, test_data.v, DisplayName='Voltage'); hold on;
    xlabel('Time (h)'); ylabel('Voltage (V)'); legend(Location="best");
    thesisFormat;
    % -------------------------------------------------------------------------
    cfg = struct('name', {'CC_dis','CV_dis','CC_chg','CV_chg'}, ... % test phases
                 'step', {5, 6, 3, 4}, ...                          % step values of phases for processing
                 'counter', {1, 1, 1, 1} ...                        % counter values of phases for processing
                 );
    voltChannels = {'MainVolt'};
    filtv = settings.voltageRange;
    [cal, ocvrel, OCVestimators] = processOCPCurves(test_data, cfg, ...
        voltChannels, 1, cellName, filtv, @(Z) 0 < Z & Z < 0.1, 2e3, temperatureDegC);
    % -------------------------------------------------------------------------
    ocpFigures(2) = figure('Name', sprintf('2. %s: OCP on relative lithiation axis', settings.electrodeName), 'NumberTitle', 'off');
    for k = 1:length(ocvrel.estimates)
        plot(1-ocvrel.estimates(k).Z, ocvrel.estimates(k).U); hold on
    end
    xlabel('Lithiation $\tilde\theta$', Interpreter='latex');
    ylabel('Potential vs. Li (V)');
    if settings.showRelativeLegend
        legend('$U_\mathrm{ocp}^\mathrm{dch}$', '$U_\mathrm{ocp}^\mathrm{chg}$', '$U_\mathrm{ocp}^\mathrm{avg}$', Location='best');
    end
    ylim(v_range); xlim(sol_range);
    % Insets use normalized figure positions and dimensionless lithiation limits.
    mainAxes = gca;
    for inset = settings.relativeInsets'
        inset_ax = createInset(mainAxes, inset);
        for k = 1:length(ocvrel.estimates)
            plot(1-ocvrel.estimates(k).Z, ocvrel.estimates(k).U, LineWidth=2); hold on
        end
        xlim(inset_ax, inset(5:6));
    end
    thesisFormat;

    %% Rough estimate of absolute lithiation
    % Preserve the original endpoint assumptions and branch choices.
    if settings.lithiationAnchor == "delithiated"
        idx = find(string({ocvrel.estimates.name}) == "charge", 1);
        ocvrel.theta_min = 0;
        ocvrel.theta_max = max(1-ocvrel.estimates(idx).Z)*ocvrel.Q/Q_r_AM;
    else
        idx = find(string({ocvrel.estimates.name}) == "discharge", 1);
        ocvrel.theta_max = 1;
        ocvrel.theta_min = 1-max(1-ocvrel.estimates(idx).Z)*ocvrel.Q/Q_r_AM;
    end
    for k = 1:length(ocvrel.estimates)                                         
        Z = 1 - ocvrel.estimates(k).Z;                                  % Change SOC to relative lithiation
        ocvrel.estimates(k).theta = ocvrel.theta_min + Z .* (ocvrel.theta_max - ocvrel.theta_min);
    end
    % -------------------------------------------------------------------------
    ocpFigures(3) = figure('Name', sprintf('3. %s: Cutoff-based absolute lithiation', settings.electrodeName), 'NumberTitle', 'off');
    for k = 1:length(ocvrel.estimates)
        plot(ocvrel.estimates(k).theta, ocvrel.estimates(k).U)
        hold on
    end
    xlabel('Lithiation $\theta$', Interpreter='latex');
    ylabel('Potential vs. Li (V)');
    legend('$U_\mathrm{ocp}^\mathrm{dch}$', '$U_\mathrm{ocp}^\mathrm{chg}$', '$U_\mathrm{ocp}^\mathrm{avg}$', Location='best'); 
    ylim(v_range); xlim(sol_range);
    mainAxes = gca;
    for inset = settings.relativeInsets'
        inset_ax = createInset(mainAxes, inset);
        for k = 1:length(ocvrel.estimates)
            plot(ocvrel.estimates(k).theta, ocvrel.estimates(k).U, LineWidth=2); hold on
        end
        % Show the same measured regions as figure 2 on the absolute axis.
        thetaLimits = ocvrel.theta_min + inset(5:6) .* ...
            (ocvrel.theta_max - ocvrel.theta_min);
        xlim(inset_ax, thetaLimits);
    end
    thesisFormat;
    drawnow
    %% MSMR model fit
    config.init = settings.initialMSMR;
    if isempty(config.init)
        [U0, X, w] = guessMSMR(ocvrel, settings.galleryCount, temperatureDegC);
        if isempty(X) || isempty(w)
            error('modelOCPelectrode:IncompleteInitialization', ...
                'Not enough peaks for the requested gallery count; inspect guessMSMR diagnostics.');
        end
        config.init = struct('U0', U0, 'X', X, 'w', w);
    end
    config.init.thetamin = 0;
    config.init.thetamax = 1;
    parameterNames = {'U0', 'X', 'w'};
    for i = 1:numel(parameterNames)
        parameter = parameterNames{i};
        config.init.(parameter) = config.init.(parameter)(:);
        config.lb.(parameter) = config.init.(parameter)*settings.boundScales(i);
        config.ub.(parameter) = config.init.(parameter)/settings.boundScales(i);
    end
    config.lb.thetamin = 0;
    config.lb.thetamax = 0.85;
    config.ub.thetamin = 0.1;
    config.ub.thetamax = 1;
    config.optw = settings.differentialWeight;
    if isfield(settings, 'optimizer')
        config.optimizer = settings.optimizer;
    end
    cons = struct('A', [], 'b', []);
    % -------------------------------------------------- Launch MSMR fit
    vis = 0;
    MSMR_curves = ["discharge", "charge", "diagonalAverage"];
    for curve = MSMR_curves
        [params, Zmod, Zrel, dZmod, xj, Uj, dZj, fun, fitInfo] = fitMSMR(config, ocvrel, temperatureDegC, curve, vis, cons); 
        ocvrel.MSMR.(curve).params = params;
        ocvrel.MSMR.(curve).fitInfo = fitInfo;
        ocvrel.MSMR.(curve).Zmod = Zmod;
        ocvrel.MSMR.(curve).Zrel = Zrel;
        ocvrel.MSMR.(curve).dZmod = dZmod;
        ocvrel.MSMR.(curve).xj = xj;
        ocvrel.MSMR.(curve).Uj = Uj;
        ocvrel.MSMR.(curve).dZj = dZj;
        ocvrel.MSMR.(curve).fun = fun;
        % Zmod was evaluated at these potentials; inversion adds no information.
        idx = find(string({ocvrel.estimates.name}) == curve, 1);
        ocvrel.MSMR.(curve).U = ocvrel.estimates(idx).U(:);
        if ~fitInfo.success
            warning('modelOCPelectrode:FitNotConverged', ...
                '%s: inspect MSMR.%s.fitInfo before interpreting the fitted calibration.', ...
                settings.electrodeName, curve);
        end
    end
    % -------------------------------------------------- Relative lithiation
    ocpFigures(4) = figure('Name', sprintf('4. %s: MSMR fits on relative lithiation axis', settings.electrodeName), 'NumberTitle', 'off');
    for curve = MSMR_curves
        plot(1-ocvrel.estimates(find(string({ocvrel.estimates.name}) == curve)).Z, ...
             ocvrel.estimates(find(string({ocvrel.estimates.name}) == curve)).U); hold on
    end
    for curve = MSMR_curves % fit
        plot(ocvrel.MSMR.(curve).Zrel, ocvrel.MSMR.(curve).U, 'k:')
    end
    ylabel('Potential vs. Li (V)');
    xlabel('Lithiation $\tilde\theta$', Interpreter='latex');
    legend('C/100 dch', 'C/100 chg', 'C/100 avg', ...
        'Model fit', '', '', Location='best')
    ylim(v_range); xlim(sol_range)
    mainAxes = gca;
    for inset = settings.fitInsets'
        inset_ax = createInset(mainAxes, inset);
        for curve = MSMR_curves
            idx = find(string({ocvrel.estimates.name}) == curve, 1);
            plot(1-ocvrel.estimates(idx).Z, ocvrel.estimates(idx).U); hold on
        end
        for curve = MSMR_curves
            plot(ocvrel.MSMR.(curve).Zrel, ocvrel.MSMR.(curve).U, 'k:')
        end
        xlim(inset_ax, inset(5:6));
    end
    thesisFormat;
    % Retain residual diagnostics without creating a separate figure.
    for curve = MSMR_curves
        idx = find(string({ocvrel.estimates.name}) == curve, 1);
        Zexp = 1 - ocvrel.estimates(idx).Z(:);
        Uexp = ocvrel.estimates(idx).U(:);
        funModel = ocvrel.MSMR.(curve).fun;
        fitParams = ocvrel.MSMR.(curve).params;
        % The inverse model expects absolute composition, not measured relative lithiation.
        thetaExp = fitParams.thetamin + ...
            (Zexp-min(Zexp))/(max(Zexp)-min(Zexp)) * ...
            (fitParams.thetamax-fitParams.thetamin);
        Umodel = zeros(size(Zexp));
        for k = 1:numel(Zexp)
            Umodel(k) = funModel(thetaExp(k));
        end
        err = Umodel - Uexp; % Voltage residual: model minus measurement.
        ocvrel.MSMR.(curve).voltageResidual = struct( ...
            'relativeLithiation', Zexp, 'voltage_V', err, ...
            'rmse_V', sqrt(mean(err.^2)), 'maxAbs_V', max(abs(err)));
    end
    % -------------------------------------------------- Average MSMR
    MSMR_curves(end+1) = "average";
    paramsChg = ocvrel.MSMR.charge.params;
    paramsDch = ocvrel.MSMR.discharge.params;
    paramsAvg.U0 = (paramsChg.U0 + paramsDch.U0) / 2;
    paramsAvg.X = (paramsChg.X + paramsDch.X) / 2;
    paramsAvg.w = (paramsChg.w + paramsDch.w) / 2;
    paramsAvg.thetamin = (paramsChg.thetamin + paramsDch.thetamin) / 2;
    paramsAvg.thetamax = (paramsChg.thetamax + paramsDch.thetamax) / 2;
    ocvrel.MSMR.average.params = paramsAvg;
    [Zmod,dZmod,xj,Uj,dZj] = compositionMSMR(ocvrel.MSMR.charge.U, paramsAvg, temperatureDegC);
    Zrel = (Zmod - paramsAvg.thetamin) ./ (paramsAvg.thetamax - paramsAvg.thetamin);
    Zrel = min(ocvrel.MSMR.charge.Zrel) + Zrel * (max(ocvrel.MSMR.charge.Zrel) - min(ocvrel.MSMR.charge.Zrel));
    U = ocvrel.MSMR.charge.U(:);
    ocvrel.MSMR.average.Zmod = Zmod;
    ocvrel.MSMR.average.Zrel = Zrel;
    ocvrel.MSMR.average.dZmod = dZmod;
    ocvrel.MSMR.average.xj = xj;
    ocvrel.MSMR.average.Uj = Uj;
    ocvrel.MSMR.average.dZj = dZj;
    ocvrel.MSMR.average.U = U;
    % -------------------------------------------------- Absolute lithiation
    ocpFigures(5) = figure('Name', sprintf('5. %s: MSMR fits on absolute lithiation axis', settings.electrodeName), 'NumberTitle', 'off');
    for curve = MSMR_curves
        plot(ocvrel.MSMR.(curve).Zmod, ocvrel.MSMR.(curve).U); hold on
    end
    ylabel('Potential vs. Li (V)');
    xlabel('Lithiation $\theta$', Interpreter='latex');
    legend('dch', 'chg', 'avg', 'mean', Location='best'); 
    ylim(v_range); xlim(sol_range);
    mainAxes = gca;
    for inset = settings.absoluteInsets'
        inset_ax = createInset(mainAxes, inset);
        for curve = MSMR_curves
            plot(ocvrel.MSMR.(curve).Zmod, ocvrel.MSMR.(curve).U); hold on
        end
        xlim(inset_ax, inset(5:6));
    end
    thesisFormat;
    drawnow

    %% MSMR model plots
    ocpFigures(6) = figure('Name', sprintf('6. %s: Average MSMR composition and galleries', settings.electrodeName), 'NumberTitle', 'off');
    subplot(1, 3, 1);
    plot(ocvrel.MSMR.average.Zmod, ocvrel.MSMR.average.U, '-',DisplayName='$U_\mathrm{ocp}^\mathrm{MSMR}$'); hold on
    plot(ocvrel.MSMR.average.dZmod./max(ocvrel.MSMR.average.dZmod), ocvrel.MSMR.average.U, '-',DisplayName='Normalized $|\mathrm{d}\theta/\mathrm{d}U|$', LineWidth=1.25);
    xlabel('Lithiation $\theta$ or normalized $|\mathrm{d}\theta/\mathrm{d}U|$', Interpreter='latex');
    ylabel('Potential vs. Li (V)');
    ylim(v_range); xlim(sol_range);
    legend(Location='best');
    % --------------------------------------------------
    subplot(1, 3, 2);
    for k = 1:size(ocvrel.MSMR.average.xj, 1)
        plot(ocvrel.MSMR.average.xj(k, :), ocvrel.MSMR.average.Uj(k, :), DisplayName=sprintf('$j$ = %d', k)); hold on
    end
    xlabel('Fractional occupation of Li $x_\mathrm{j}$', Interpreter='latex');
    ylabel('Potential vs. Li (V)');
    ylim(v_range);
    legend(Location='best');
    % --------------------------------------------------
    subplot(1, 3, 3);
    for k = 1:size(ocvrel.MSMR.average.dZj, 1)
        plot(ocvrel.MSMR.average.Uj(k, :), thermalFactor*ocvrel.MSMR.average.dZj(k, :), DisplayName=sprintf('$j$ = %d', k)); hold on
    end
    xlabel('Potential vs. Li (V)');
    ylabel('$|\mathrm{d}x_j/\mathrm{d}U|$ (V$^{-1}$)', Interpreter='latex');
    xlim(v_range);
    legend(Location='best');
    thesisFormat;
    drawnow

    %% Supplementary S3.3: retain results before AM/hidden-gallery corrections
    % Value copies keep the fitted intervals and galleries available after later edits.
    supplementaryResults = struct('electrode', settings.electrodeName, ...
        'temperatureDegC', temperatureDegC, 'raw', ocvrel.estimates, 'MSMR', ocvrel.MSMR, ...
        'configuration', config, 'halfCellCharge_Ah', ocvrel.Q, ...
        'initialCapacity_Ah', Q_r_AM, 'figures', ocpFigures(1:6));
    %% Optional hidden-gallery correction (after the supplementary snapshot)
    % Retain the original transformation; the snapshot above precedes it.
    if settings.hiddenGalleryCorrection
        for curve = MSMR_curves
            ocvrel.MSMR.(curve).params.thetamin = ocvrel.theta_min;
            ocvrel.MSMR.(curve).Zmod = ocvrel.MSMR.(curve).Zmod .* ...
                (ocvrel.MSMR.(curve).params.thetamax - ocvrel.MSMR.(curve).params.thetamin) + ...
                ocvrel.MSMR.(curve).params.thetamin;
        end
    end

    %% AM correction (after the supplementary snapshot)
    Q_MSMR      = ocvrel.Q/(ocvrel.MSMR.average.params.thetamax - ocvrel.MSMR.average.params.thetamin);
    AM_per_diff = Q_MSMR/Q_r - AM_per;
    AM_per      = AM_per + AM_per_diff;
    Q_r_AM      = Q_r*AM_per;                          % corrected sample capacity
    % --------------------------------------------------
    idx = find(string({ocvrel.estimates.name}) == "diagonalAverage");
    ocpFigures(7) = figure('Name', sprintf('7. %s: Relative, cutoff-based, and MSMR average OCP', settings.electrodeName), 'NumberTitle', 'off');
    plot(1-ocvrel.estimates(idx).Z, ocvrel.estimates(idx).U, DisplayName='$U_\mathrm{ocp}\bigl(\tilde\theta^\mathrm{r}\bigr)$'); hold on
    plot(ocvrel.estimates(idx).theta, ocvrel.estimates(idx).U, DisplayName='$U_\mathrm{ocp}(\theta^\mathrm{r})$');
    plot(ocvrel.MSMR.average.Zmod, ocvrel.MSMR.average.U, DisplayName='$U_\mathrm{ocp}^\mathrm{MSMR}$');
    xlabel('Lithiation $\tilde\theta,\ \theta$', Interpreter='latex');
    ylabel('Potential vs. Li (V)');
    ylim(v_range); xlim(sol_range);
    if settings.showComparisonLegend
        legend(Location='best');
    end
    mainAxes = gca;
    for inset = settings.relativeInsets'
        inset_ax = createInset(mainAxes, inset);
        plot(1-ocvrel.estimates(idx).Z, ocvrel.estimates(idx).U, LineWidth=2); hold on
        plot(ocvrel.estimates(idx).theta, ocvrel.estimates(idx).U, LineWidth=2);
        plot(ocvrel.MSMR.average.Zmod, ocvrel.MSMR.average.U, LineWidth=2);
        % Reuse figure 3's inset regions on the absolute lithiation axis.
        thetaLimits = ocvrel.theta_min + inset(5:6) .* ...
            (ocvrel.theta_max - ocvrel.theta_min);
        xlim(inset_ax, thetaLimits);
    end
    thesisFormat;
    drawnow
    % --------------------------------------------------
    for k = 1:length(ocvrel.estimates)  % Change absolut lithiation range to match MSMR
        ocvrel.estimates(k).theta = ocvrel.MSMR.average.params.thetamin + (1 - ocvrel.estimates(k).Z) .* (ocvrel.MSMR.average.params.thetamax - ocvrel.MSMR.average.params.thetamin);
    end

    %% Assemble the reusable model; save only when explicitly requested
    % ---- Uocp estimates
    model.raw = ocvrel.estimates;
    % ---- MSMR
    model.MSMR = ocvrel.MSMR;
    % ---- LUT
    model.LUT.theta_ch     = model.raw(find(string(OCVestimators) == "charge")).theta;
    model.LUT.thetarel_ch  = (1 - model.raw(find(string(OCVestimators) == "charge")).Z);
    model.LUT.U_ch         = model.raw(find(string(OCVestimators) == "charge")).U;
    model.LUT.UocpCh       = @(theta) linearinterp(model.LUT.theta_ch, model.LUT.U_ch, theta);

    model.LUT.theta_dch    = model.raw(find(string(OCVestimators) == "discharge")).theta;
    model.LUT.thetarel_dch = (1 - model.raw(find(string(OCVestimators) == "discharge")).Z);
    model.LUT.U_dch        = model.raw(find(string(OCVestimators) == "discharge")).U;
    model.LUT.UocpDch      = @(theta) linearinterp(model.LUT.theta_dch, model.LUT.U_dch, theta);

    model.LUT.theta_avg    = model.raw(find(string(OCVestimators) == "diagonalAverage")).theta;
    model.LUT.thetarel_avg = (1 - model.raw(find(string(OCVestimators) == "diagonalAverage")).Z);
    model.LUT.U_avg        = model.raw(find(string(OCVestimators) == "diagonalAverage")).U;
    model.LUT.UocpAvg      = @(theta) linearinterp(model.LUT.theta_avg, model.LUT.U_avg, theta);
    % ---- fcn
    model.zToTheta = @(z) z .* (model.MSMR.average.params.thetamin - model.MSMR.average.params.thetamax) + model.MSMR.average.params.thetamax;
    model.thetaToZ = @(theta) (theta - model.MSMR.average.params.thetamax) ./ (model.MSMR.average.params.thetamin - model.MSMR.average.params.thetamax);
    % --------------------------------------------------
    model.Q25 = ocvrel.Q;
    model.Qnom = Q_r_AM;
    % Retain the experimental assumptions and fit settings alongside the LUT.
    model.provenance = struct('settings', settings, 'configuration', config, ...
        'matlabVersion', version, 'calibration', 'MSMR interval after configured corrections');
    % --------------------------------------------------
    if strlength(string(settings.outputFile)) > 0
        save(settings.outputFile, 'model');
    end
end

function insetAxes = createInset(mainAxes, inset)
%CREATEINSET Create an inset associated with its main axes for formatting.
    inset = inset(:).';
    insetAxes = axes('Parent', ancestor(mainAxes, 'figure'), ...
        'Units', 'normalized', 'Position', inset(1:4), 'Box', 'on');
    setappdata(insetAxes, 'thesisFormatMainAxes', mainAxes);
end
