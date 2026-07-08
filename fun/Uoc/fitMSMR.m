function [paramsfin, Zmodfin, Zrelfin, dZmodfin, xj, Uj, dZj, fun] = fitMSMR(config, ocvrel, TdegC, curve, vis, cons)

    J      = length(config.init.U0); % number of MSMR galleries
    indU   = 1:J;                    % indicies to U0    in param vectors
    indX   = J+1:2*J;                % indicies to X     in param vectors
    indW   = 2*J+1:3*J;              % indicies to omega in param vectors
    % indX = J+1:2*J-1;              % indicies to X     in param vectors
    % indW = 2*J:3*J-1;              % indicies to omega in param vectors
    weight = config.optw;            % Relative weight to assign to dZdV residuals
    Usep   = 0.001;                  % Minimum voltage separation between peaks (1mV)
    
    % ------------------------------------------------------- Find the OCP vector to fit MSMR
    idx      = find(string({ocvrel.estimates.name}) == curve);                               
    Ulab     = ocvrel.estimates(idx).U;                        % OCP vector derived from lab data [V]
    Zlab     = 1-ocvrel.estimates(idx).Z;                      % SOC vector derived from lab data [-] (Relative lithiation)
    dZlab    = abs(ocvrel.estimates(idx).dZ);                  % Absolute differential capacity
    diffzlab = max(Zlab) - min(Zlab);                          % Relative lithiation range

    % -------------------------------------------------------
    [Zmod, dZmod] = compositionMSMR(Ulab, config.init, TdegC);                                   % Evaluate MSMR
    Zmod          = (Zmod - config.init.thetamin)/(config.init.thetamax - config.init.thetamin); % From absolute to relative
    Zmod          = min(Zlab) + Zmod*diffzlab;
    dZmod         = dZmod*diffzlab/(config.init.thetamax - config.init.thetamin);                % From absolute to relative
    
    % -------------------------------------------------------
    if vis == 1
        persistent fig;
        fig = figure; clf(fig);
        subplot(2,2,1)
            line_OcpLab = plot(Zlab,Ulab,'-', DisplayName='Lab'); hold on;
            line_OcpIni = plot(Zmod, Ulab, 'r:', DisplayName='Ini');
            line_OcpMod = plot(NaN,NaN,'k:', DisplayName='Model');
            xlim([0 1])
            ylabel('$U_\mathrm{ocp}$ ($\mathrm{V}$)','Interpreter','latex');
            legend(Location="best")
        subplot(2,2,2)
            line_dOcpLab = plot(dZlab./max(dZlab),Ulab,'-', DisplayName='Lab'); hold on;
            line_dOcpIni = plot(dZmod./max(dZmod), Ulab, 'r:', DisplayName='Ini');
            line_dOcpMod = plot(NaN,NaN,'k:', DisplayName='Model');
            ylim([min(config.init.U0)-0.25 max(config.init.U0)+0.25])
            xlabel('Normalized $d\tilde{\theta}_\mathrm{s}/dU$', 'Interpreter','latex');
        subplot(2,2,3)
            line_dthetaLab = plot(Zlab,dZlab./max(dZlab),'-', DisplayName='Lab'); hold on;
            line_dthetaIni = plot(Zmod, dZmod./max(dZmod), 'r:', DisplayName='Ini');
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
        % params.X = params.X(1:end-1);
        x = [params.U0; params.X; params.w; params.thetamin; params.thetamax];
    end

    function params = vectToParams(x, indU, indX, indW)
        x = x(:);
        params.U0 = x(indU);  
        % params.X = [x(indX); 1-sum(x(indX))]; % Same to applying cons
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
    orderU.U0       = toeplitz([-1 zeros(1,J-2)],[-1 1 zeros(1,J-2)])';
    orderU.X        = zeros(J,J-1);
    orderU.w        = zeros(J,J-1);
    orderU.thetamin = zeros(1,J-1);
    orderU.thetamax = zeros(1,J-1);
    A               = -paramsToVect(orderU)';     % U1 < U2 < U3...
    b               = -ones(J-1,1).*Usep(:);      % minimum voltage peak separation of Usep;
    A               = [A; cons.A];                % Addition of additional cons
    b               = [b; cons.b];

    % ------------------------------------------------------- Non-linear cons.
    function [c, ceq] = nonlcon(x,indU, indX, indW, TdegC, Zlab, diffzlab, Umod) % Conversion to relative must be 0-1 range
        params = vectToParams(x, indU, indX, indW);                                 % Unpack parameter vector
        [Zmod, ~] = compositionMSMR(Umod, params, TdegC);                           % Evaluate MSMR  
        Zmod = (Zmod - params.thetamin)/(params.thetamax - params.thetamin);        % From absolute to relative, Eq.: 3.2
        Zmod = min(Zlab) + Zmod*diffzlab;
        c = (max(Zmod) - min(Zmod)) - 1;
        ceq = [];
    end

    % ------------------------------------------------------- Optimizer opts.
    optsHybrid = optimoptions(@fmincon, ...
      'Display','iter',...
      'MaxFunEvals',1e6,...   % default: 3000
      'MaxIter',1000,...      % default: 1000
      'TolFun',1e-10,...      % default: 1e-6
      'TolX',1e-10,...        % default: 1e-10
      'TolCon',1e-10...       % default: 1e-6
    );
    opts = optimoptions(@ga,...
        'Display', 'iter', ...
        'PopulationSize', 200, ...                  % default: 200
        'MaxGenerations', 20, ...                   % default: 100*numberOfVariables
        'FunctionTolerance', 1e-10,...              % default: 1e-6
        'MaxStallGenerations', 100,...              % default: 50
        'MaxTime', 120,...
        'HybridFcn', {@fmincon,optsHybrid}, ...
        'InitialPopulationMatrix', paramsToVect(config.init)', ...
        'UseParallel', true ...
    );

    % ------------------------------------------------------- Optimization
    fun = @(x) fitMSMR_costFcn(x, Ulab, TdegC, Zlab, diffzlab, dZlab, weight, indU, indX, indW, vis);
    % x   = ga(fun,...
    %          length(paramsToVect(config.init)), ...
    %          A,b,Aeq,beq,lb,ub, ...
    %          @(x) nonlcon(x, indU, indX, indW, TdegC, Zlab, diffzlab, Ulab), ...
    %          opts); % Hybrid optimization
    x = fmincon(fun,paramsToVect(config.init),A,b,Aeq,beq,lb,ub,@(x) nonlcon(x, indU, indX, indW, TdegC, Zlab, diffzlab, Ulab),optsHybrid);
    
    % ------------------------------------------------------- Model
    paramsfin = vectToParams(x, indU, indX, indW);                                        % Unpack parameters to vector
    [Zmodfin, dZmodfin, xj, Uj, dZj] = compositionMSMR(Ulab, paramsfin, TdegC);           % Evaluate MSMR
    Zrelfin = (Zmodfin - paramsfin.thetamin)/(paramsfin.thetamax - paramsfin.thetamin);   % From absolute to relative, Eq.: 3.2
    Zrelfin = min(Zlab) + Zrelfin*diffzlab;
    fun = @(z) fzero(@(U) lithiationBalanceEquation(U, z, paramsfin, TdegC), mean(Ulab)); % solving the nonlinear equation in Uocp

    %% ------------------------------------------------------- Cost fcn.
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