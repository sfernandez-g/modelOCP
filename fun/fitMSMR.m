% FITMSMR Fit a multi-site multi-reaction (MSMR) model to experimental OCP data.
%
% -- Usage --
% [paramsfin,Zmodfin,Zrelfin,dZmodfin,xj,Uj,dZj,fun,fitInfo] = ...
%     fitMSMR(config,ocvrel,TdegC,curve,vis,cons)
%
% -- Description --
% Fits an MSMR model to an experimental open-circuit potential (OCP) curve.
% The model parameters are obtained using a two-stage optimization:
%
%   1. Global search using a genetic algorithm (GA).
%   2. Local refinement of the GA solution using fmincon.
%
% The MSMR model is evaluated as composition Z(U), where each gallery is
% defined by its equilibrium potential U0, capacity fraction X, and
% transition-width parameter w. The fractions X are constrained to sum to 1,
% and consecutive U0 values are constrained to remain separated by a minimum
% voltage interval.
%
% The experimental OCP curve is selected from ocvrel.estimates using "curve".
% The experimental composition is converted to relative lithiation before
% fitting. The model composition is transformed using the fitted
% thetamin and thetamax parameters to match the experimental composition
% range.
%
% -- Inputs --
% config   = Structure containing:
%            .init.U0       Initial values of MSMR equilibrium potentials [V].
%            .init.X        Initial gallery capacity fractions [-].
%            .init.w        Initial transition-width parameters.
%            .init.thetamin Initial minimum composition.
%            .init.thetamax Initial maximum composition.
%            .lb            Lower parameter bounds.
%            .ub            Upper parameter bounds.
%            .optw          Weight applied to the differential-composition
%                           residual.
%
% ocvrel  = Structure containing experimental/estimated OCP curves.
%
% TdegC   = Temperature of the experiment [degC].
%
% curve   = Name of the OCP curve in ocvrel.estimates to be fitted.
%
% vis     = Logical flag. If 1, display the fitting process and parameter
%           evolution during optimization.
%
% cons    = Structure containing additional linear inequality constraints:
%            .A
%            .b
%           These are appended to the MSMR U0-ordering constraints.
%
% -- Outputs --
% paramsfin = Fitted MSMR parameters:
%             .U0
%             .X
%             .w
%             .thetamin
%             .thetamax
%
% Zmodfin  = MSMR composition evaluated at the experimental voltage points.
%
% Zrelfin  = MSMR composition transformed to the relative experimental
%             lithiation range.
%
% dZmodfin = Differential MSMR composition dZ/dU.
%
% xj       = Composition contribution of each MSMR gallery.
%
% Uj       = Equilibrium potential associated with each gallery contribution.
%
% dZj      = Differential contribution of each MSMR gallery.
%
% fun      = Function handle that calculates the OCP voltage corresponding
%            to a specified relative composition.
%
% -- Optimization --
% The objective function minimizes the combined squared residuals between
% experimental and model composition and differential composition:
%
%   J = sum((Zmod-Zlab).^2) + optw*sum((dZmod-dZlab).^2)
%
% The optimization is subject to parameter bounds, sum(X)=1, the ordering
% and minimum separation of consecutive U0 values, and any additional
% constraints supplied through cons.

function [paramsfin, Zmodfin, Zrelfin, dZmodfin, xj, Uj, dZj, fun, fitInfo] = fitMSMR(config, ocvrel, TdegC, curve, vis, cons)

% Optional optimizer settings: config.optimizer.seed, .gaOptions, and
% .fminconOptions. The options fields accept name/value cell arrays;
% omitted settings retain the defaults below.
optimizer = struct('seed', 0, 'gaOptions', {{}}, 'fminconOptions', {{}});
if isfield(config, 'optimizer')
    supplied = config.optimizer;
    fields = fieldnames(supplied);
    for k = 1:numel(fields)
        if isfield(optimizer, fields{k})
            optimizer.(fields{k}) = supplied.(fields{k});
        end
    end
end
if ~isscalar(optimizer.seed) || ~isfinite(optimizer.seed) || optimizer.seed < 0
    error('fitMSMR:InvalidSeed', 'config.optimizer.seed must be a finite nonnegative scalar.');
end
rngState = rng;
rngCleanup = onCleanup(@() rng(rngState)); %#ok<NASGU>
rng(optimizer.seed, 'twister');

    J      = length(config.init.U0); % number of MSMR galleries
    indU   = 1:J;                    % indicies to U0    in param vectors
    indX   = J+1:2*J;                % indicies to X     in param vectors
    indW   = 2*J+1:3*J;              % indicies to omega in param vectors
    weight = config.optw;            % Relative weight to assign to dZdV residuals
    Usep   = 0.001;                  % Minimum voltage separation between peaks (1mV)
    
    % ------------------------------------------------------- Find the OCP vector to fit MSMR
    idx      = find(string({ocvrel.estimates.name}) == curve);                               
    Ulab     = ocvrel.estimates(idx).U;                        % OCP vector derived from lab data [V]
    Zlab     = 1-ocvrel.estimates(idx).Z;                      % SOC vector derived from lab data [-] (Relative lithiation)
    dZlab    = abs(ocvrel.estimates(idx).dZ);                  % Absolute differential capacity
    diffzlab = max(Zlab) - min(Zlab);                          % Relative lithiation range

    % -------------------------------------------------------
    [Zmod_ini, dZmod_ini] = compositionMSMR(Ulab, config.init, TdegC);                                   % Evaluate MSMR
    Zmod_ini          = (Zmod_ini - config.init.thetamin)/(config.init.thetamax - config.init.thetamin); % From absolute to relative
    Zmod_ini          = min(Zlab) + Zmod_ini*diffzlab;
    dZmod_ini         = dZmod_ini*diffzlab/(config.init.thetamax - config.init.thetamin);                % From absolute to relative
    
    % -------------------------------------------------------
    if vis == 1
        persistent fig;
        fig = figure; clf(fig);
        subplot(2,2,1)
            line_OcpLab = plot(Zlab,Ulab,'-', DisplayName='Lab'); hold on;
            line_OcpIni = plot(Zmod_ini, Ulab, 'r:', DisplayName='Ini');
            line_OcpMod = plot(NaN,NaN,'k:', DisplayName='Model');
            xlim([0 1])
            ylabel('$U_\mathrm{ocp}$ ($\mathrm{V}$)','Interpreter','latex');
            legend(Location="best")
        subplot(2,2,2)
            line_dOcpLab = plot(dZlab./max(dZlab),Ulab,'-', DisplayName='Lab'); hold on;
            line_dOcpIni = plot(dZmod_ini./max(dZmod_ini), Ulab, 'r:', DisplayName='Ini');
            line_dOcpMod = plot(NaN,NaN,'k:', DisplayName='Model');
            ylim([min(config.init.U0)-0.25 max(config.init.U0)+0.25])
            xlabel('Normalized $d\tilde{\theta}_\mathrm{s}/dU$', 'Interpreter','latex');
        subplot(2,2,3)
            line_dthetaLab = plot(Zlab,dZlab./max(dZlab),'-', DisplayName='Lab'); hold on;
            line_dthetaIni = plot(Zmod_ini, dZmod_ini./max(dZmod_ini), 'r:', DisplayName='Ini');
            line_dthetaMod = plot(NaN,NaN,'k:', DisplayName='Model');
            xlabel('Relative Composition, $\tilde{\theta}_\mathrm{s}$', 'Interpreter','latex');
            ylabel('Normalized $d\tilde{\theta}_\mathrm{s}/dU$', 'Interpreter','latex');
            xlim([0 1])
        
        ax4 = subplot(2,2,4); pos = get(ax4, 'Position'); delete(ax4)
        for i = 1:3
            ax(i) = axes('Position', [pos(1), pos(2) + (3-i)*(pos(4)/3) + 0.05/i, ...
                                     pos(3), ...
                                     pos(4)/3/1.5]);
        end   
        axes(ax(1))
            line_U0 = plot(indU,config.init.U0,'r:', 'Marker','o', 'DisplayName','$U_\mathrm{j}^\mathrm{0}$'); hold on;
            line_U0 = plot(NaN,NaN,'k:', 'Marker','o', 'DisplayName','$U_\mathrm{j}^\mathrm{0}$');
            ylabel('$U_j^0$', 'Interpreter','latex')   
        axes(ax(2))
            line_X = plot(indU,config.init.X,'r:', 'Marker','o', 'DisplayName','$X_\mathrm{j}$'); hold on;
            line_X = plot(NaN,NaN,'k:', 'Marker','o', 'DisplayName','$X_\mathrm{j}$');
            ylabel('$X_j$', 'Interpreter','latex')
        axes(ax(3))
            line_w = plot(indU,config.init.w,'r:', 'Marker','o', 'DisplayName','$w_\mathrm{j}$'); hold on;
            line_w = plot(NaN,NaN,'k:', 'Marker','o', 'DisplayName','$w_\mathrm{j}$');
            ylabel('$w_j$', 'Interpreter','latex')
            xlabel('Gallery $j$', 'Interpreter','latex')
        for i = 1:3
            xlim(ax(i), [min(indU) max(indU)])
        end
        set(ax(1), 'XTickLabel', [])
        set(ax(2), 'XTickLabel', [])
    end
    
    % ------------------------------------------------------- Array flips
    function x = paramsToVect(params)
        x = [params.U0; params.X; params.w; params.thetamin; params.thetamax];
    end
    function params = vectToParams(x, indU, indX, indW)
        x = x(:);
        params.U0 = x(indU);  
        params.X = x(indX); 
        params.w = x(indW);
        params.thetamin = x(end-1); 
        params.thetamax = x(end); 
    end

    % ------------------------------------------------------- Bounds
    lb = paramsToVect(config.lb);   
    ub = paramsToVect(config.ub);

    % ------------------------------------------------------- Equality cons.
    Aeq       = zeros(1, length(lb)); 
    Aeq(indX) = 1; 
    beq       = 1;  % sum of Xj equal to 1
    
    % ------------------------------------------------------- Non-equality cons.
    orderU = zeros(J-1, length(lb));
    for j = 1:J-1
        orderU(j, indU(j))   = 1;
        orderU(j, indU(j+1)) = -1;
    end
    A = orderU;                   % U0(j+1) - U0(j) >= Usep
    b = -Usep * ones(J-1, 1);
    A = [A; cons.A];
    b = [b; cons.b];

    % ------------------------------------------------------- Optimizer opts.
    optsGA = optimoptions(@ga, ...
        'Display', 'iter', ...
        'PopulationSize', 100, ...
        'MaxGenerations', 100, ...
        'MaxStallGenerations', 5, ...
        'FunctionTolerance', 1e-8, ...
        'ConstraintTolerance', 1e-8, ...
        'InitialPopulationMatrix',paramsToVect(config.init)',...
        'UseParallel', false);
    if ~isempty(optimizer.gaOptions)
        optsGA = optimoptions(optsGA, optimizer.gaOptions{:});
    end
    
    optsFmincon = optimoptions(@fmincon, ...
        'Display', 'iter', ...
        'MaxFunctionEvaluations', 1e6, ...
        'MaxIterations', 1000, ...
        'OptimalityTolerance', 1e-10, ...
        'StepTolerance', 1e-10, ...
        'ConstraintTolerance', 1e-10);
    if ~isempty(optimizer.fminconOptions)
        optsFmincon = optimoptions(optsFmincon, optimizer.fminconOptions{:});
    end

    % ------------------------------------------------------- Cost fcn.
    fun = @(x) fitMSMR_costFcn(x, Ulab, TdegC, Zlab, diffzlab, dZlab, weight, indU, indX, indW, 0);

    % ------------------------------------------------------- Optimization
    % Global search with GA
    [xGA, fGA, gaExitflag, gaOutput] = ga( ...
        fun, ...
        length(lb), ...
        A, b, ...
        Aeq, beq, ...
        lb, ub, ...
        [], ...
        optsGA);
    % Local refinement with fmincon
    fun = @(x) fitMSMR_costFcn(x, Ulab, TdegC, Zlab, diffzlab, dZlab, weight, indU, indX, indW, vis);
    [x, fval, fminconExitflag, fminconOutput] = fmincon( ...
        fun, ...
        xGA, ...
        A, b, ...
        Aeq, beq, ...
        lb, ub, ...
        [], ...
        optsFmincon);

    fitInfo.seed = optimizer.seed;
    fitInfo.ga.objective = fGA;
    fitInfo.ga.exitflag = gaExitflag;
    fitInfo.ga.output = gaOutput;
    fitInfo.fmincon.objective = fval;
    fitInfo.fmincon.exitflag = fminconExitflag;
    fitInfo.fmincon.output = fminconOutput;
    fitInfo.success = gaExitflag > 0 && fminconExitflag > 0;

    % ------------------------------------------------------- Model
    paramsfin = vectToParams(x, indU, indX, indW);                                        % Unpack parameters to vector
    [Zmodfin, dZmodfin, xj, Uj, dZj] = compositionMSMR(Ulab, paramsfin, TdegC);           % Evaluate MSMR
    Zrelfin = (Zmodfin - paramsfin.thetamin)/(paramsfin.thetamax - paramsfin.thetamin);   % From absolute to relative, Eq.: 3.2
    Zrelfin = min(Zlab) + Zrelfin*diffzlab;
    fun = @(z) fzero(@(U) lithiationBalanceEquation(U, z, paramsfin, TdegC), mean(Ulab)); % solving the nonlinear equation in Uocp

    function y = fitMSMR_costFcn(x, Umod, TdegC, Zlab, diffzlab, dZlab, weight, indU, indX, indW, vis)
        params        = vectToParams(x, indU, indX, indW);                            % Unpack parameter vector
        [Zmod, dZmod] = compositionMSMR(Umod, params, TdegC);                         % Evaluate MSMR  
        Zmod          = (Zmod - params.thetamin)/(params.thetamax - params.thetamin); % From absolute to relative
        Zmod          = min(Zlab) + Zmod*diffzlab;
        dZmod         = dZmod*diffzlab/(params.thetamax - params.thetamin);           % From absolute to relative
        
        zResid = Zmod(:) - Zlab(:);     % Compute residuals
        dzResid = dZmod(:) - dZlab(:);
        
        if vis == 1
            figure(fig)
            set(line_OcpLab,'XData',Zlab,'Ydata',Umod);
            set(line_OcpMod,'XData',Zmod,'Ydata',Umod);
            set(line_dOcpLab,'XData',dZlab./max(dZlab),'Ydata',Umod);
            set(line_dOcpMod,'XData',dZmod./max(dZmod),'Ydata',Umod);
            set(line_dthetaLab,'XData',Zlab,'Ydata',dZlab./max(dZlab));
            set(line_dthetaMod,'XData',Zmod,'Ydata',dZmod./max(dZmod));
            set(line_U0,'XData',indU,'Ydata',params.U0);
            set(line_X,'XData',indU,'Ydata',params.X);
            set(line_w,'XData',indU,'Ydata',params.w);
        end
        y = sum(zResid.^2) + weight*sum(dzResid.^2);
    end
end
