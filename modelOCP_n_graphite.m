clc; clear; close all; 
RGB = orderedcolors("gem");

%% add utility paths
addpath("fun\");

%% process OCP data
file      = 'input\neg\checkUp_LGMH1_PAT20261302_G1_Li_T25_Channel_2_Wb_1.xlsx';
%--------------------------------------------------------------------------
cellName  = 'LG MH1 PAT half-cell NE'; 
%--------------------------------------------------------------------------
v_range   = [0 1.5];
% -------------------------------------------------------------------------
sol_range   = [0 1];
% -------------------------------------------------------------------------
load("output\theoPracCapNeg.mat");
Q_r     = Q_r*1e-3;
AM_per  = 0.98;          % Rojas 2024
Q_r_AM  = Q_r*AM_per;
% -------------------------------------------------------------------------
test_data = getDataArbin_XLSX(file, 7); % filename, step where TC increments
% -------------------------------------------------------------------------
figure;
plot(test_data.t./3600, test_data.v, DisplayName='Voltage'); hold on;
pptFormatMultiplePlots('Times New Roman')
xlabel('Time (h)'); ylabel('Voltage (V)'); legend(Location="best");
% -------------------------------------------------------------------------
cfg = struct('name', {'CC_dis','CV_dis','CC_chg','CV_chg'}, ... % test phases
             'step', {5, 6, 3, 4}, ...                          % step values of phases for processing
             'counter', {1, 1, 1, 1} ...                        % counter values of phases for processing
             );
voltChannels = {'MainVolt'};
filtv        = [0, 1.5];
cal = processCU(test_data, cfg, voltChannels, 1); % test data, processing config, voltChannels

%% Estimate OCP
cellData.raw = cal; cellData.cal = cal; cellData.name = cellName;
% -------------------------------------------------------------------------
config.vmin         = cal.minV; config.vmax = cal.maxV; config.datype = 'useAvg';
config.daxcorrfiltz = @(Z){0<Z & Z<0.1};                                            % focus on this SOC region
config.daxcorrfiltv = @(V){filtv(1, 1)<V & V<filtv(1, 2)};                          % focus on this V region
% -------------------------------------------------------------------------
OCVresolution = 1e4; % this many points between 0 and 1
OCVestimators = {'diagonalAverage'};
% -------------------------------------------------------------------------
ocvrel        = estimateOCV(cellData, config, OCVresolution, OCVestimators);
ocvrel        = addDchChOCV(ocvrel, cellData, "discharge");
ocvrel        = addDchChOCV(ocvrel, cellData, "charge");
ocvrel.Q      = cal.QAh;
% -------------------------------------------------------------------------
OCVestimators{end+1} = 'charge';
OCVestimators{end+1} = 'discharge';
ocvrel.estimates     = ocvrel.estimates(end:-1:1);
% -------------------------------------------------------------------------
figure;
for k = 1:length(ocvrel.estimates)
    plot(1-ocvrel.estimates(k).Z, ocvrel.estimates(k).U); hold on
end
pptFormatMultiplePlots('Times New Roman')
xlabel('Relative degree of lithiation $\tilde{\theta}_\mathrm{s}$', Interpreter='latex'); 
ylabel('Potential vs. Li (V)');
legend('$U_\mathrm{ocp}^\mathrm{chg}$', '$U_\mathrm{ocp}^\mathrm{dch}$', '$U_\mathrm{ocp}^\mathrm{avg}$', Location='best'); 
legend(Interpreter="latex", FontName='Times New Roman', FontSize=16, Box='off', Location='best', IconColumnWidth=15)
ylim(v_range); xlim(sol_range);
title('OCP vs. Relative degree of lithiation');
drawnow
% -------------------------------------------------------------------------
figure;
for k = 1:length(ocvrel.estimates)
    plot(ocvrel.estimates(k).U, ocvrel.estimates(k).dZ); hold on
end
pptFormatMultiplePlots('Times New Roman')
ylabel('d$\tilde{\theta}_\mathrm{s}$/dU', Interpreter='latex'); 
xlabel('Potential vs. Li (V)');
legend('$U_\mathrm{ocp}^\mathrm{chg}$', '$U_\mathrm{ocp}^\mathrm{dch}$', '$U_\mathrm{ocp}^\mathrm{avg}$', Location='best'); 
legend(Interpreter="latex", FontName='Times New Roman', FontSize=16, Box='off', Location='best', IconColumnWidth=15)
xlim(v_range);
title('Differential SOL vs OCP');
drawnow
% -------------------------------------------------------------------------
figure;
for k = 1:length(ocvrel.estimates)
    plot(1-ocvrel.estimates(k).Z, ocvrel.estimates(k).dZ); hold on
end
pptFormatMultiplePlots('Times New Roman')
ylabel('d$\tilde{\theta}_\mathrm{s}$/dU', Interpreter='latex'); 
xlabel('Relative degree of lithiation $\tilde{\theta}_\mathrm{s}$', Interpreter='latex');
legend('$U_\mathrm{ocp}^\mathrm{chg}$', '$U_\mathrm{ocp}^\mathrm{dch}$', '$U_\mathrm{ocp}^\mathrm{avg}$', Location='best'); 
legend(Interpreter="latex", FontName='Times New Roman', FontSize=16, Box='off', Location='best', IconColumnWidth=15)
xlim(sol_range);
title('Differential SOL vs SOL');
drawnow

%% Rough estimate of absolute lithiation
ocvrel.theta_max = max((1-ocvrel.estimates(2).Z)).*ocvrel.Q/Q_r_AM; % Based on "discharge" data asumming it starts from a completely delithiated state
ocvrel.theta_min = 0;                                               % Assume we reach minimum lithiation
for k = 1:length(ocvrel.estimates)                                         
    Z = 1 - ocvrel.estimates(k).Z;                                  % Change SOC to relative lithiation
    ocvrel.estimates(k).theta = ocvrel.theta_min + Z .* (ocvrel.theta_max - ocvrel.theta_min);
end
% -------------------------------------------------------------------------
figure;
for k = 1:length(ocvrel.estimates)
    plot(ocvrel.estimates(k).theta, ocvrel.estimates(k).U)
    hold on
end
pptFormatMultiplePlots('Times New Roman')
xlabel('Absolute degree of lithiation $\theta_\mathrm{s}$', Interpreter='latex'); 
ylabel('Potential vs. Li (V)');
legend('$U_\mathrm{ocp}^\mathrm{chg}$', '$U_\mathrm{ocp}^\mathrm{dch}$', '$U_\mathrm{ocp}^\mathrm{avg}$', Location='best'); 
legend(Interpreter="latex", FontName='Times New Roman', FontSize=16, Box='off', Location='best', IconColumnWidth=15)
ylim(v_range); xlim(sol_range);
title('OCP vs. Absolute degree of lithiation');
drawnow

%% MSMR model fit
% --------------------------------------------------- MSMR Guess
% [U0, X, w, diag] = guessMSMR(ocvrel, 6, 25);
% config.init.U0 = U0;
% config.init.X = X;
% config.init.w = w;
% --------------------------------------------------- Literature values
config.init.U0 = [0.08843;0.12799;0.14331;0.16984;0.21446;0.36325]; % Verbrugge 2017
config.init.X = [0.43336;0.23963;0.15018;0.05462;0.06744;0.0548];
config.init.w = [0.08611;0.08009;0.72469;2.53277;0.09470;5.97354];
config.init.thetamin = 0;
config.init.thetamax = 1;
% ---------------------------------------------------
len = length(config.init.U0);
% -------------------------------------------------- Literature based bounds
config.lb.U0 = config.init.U0.*0.5;
config.lb.X = config.init.X.*0.4;
config.lb.w = config.init.w.*0.2;
config.lb.thetamin = 0-1e-10;
config.lb.thetamax = 0.85;

config.ub.U0 = config.init.U0./0.5;
config.ub.X = config.init.X./0.4;
config.ub.w = config.init.w./0.2;
config.ub.thetamin = 0.1;
config.ub.thetamax = 1+1e-10;
% -------------------------------------------------- Aditional bounds or constraints
A = []; % n x nvars
b = [];  % n x 1
cons.A = A;
cons.b = b;
% -------------------------------------------------- Differential curve weight
config.optw = 1e-3;
% -------------------------------------------------- Launch MSMR fit
vis = 1;
MSMR_curves = ["discharge", "charge", "diagonalAverage"];
for curve = MSMR_curves
    [params, Zmod, Zrel, dZmod, xj, Uj, dZj, fun] = fitMSMR(config, ocvrel, 25, curve, vis, cons); 
    ocvrel.MSMR.(curve).params = params;
    ocvrel.MSMR.(curve).Zmod = Zmod;
    ocvrel.MSMR.(curve).Zrel = Zrel;
    ocvrel.MSMR.(curve).dZmod = dZmod;
    ocvrel.MSMR.(curve).xj = xj;
    ocvrel.MSMR.(curve).Uj = Uj;
    ocvrel.MSMR.(curve).dZj = dZj;
    ocvrel.MSMR.(curve).fun = fun;
    ocvrel.MSMR.(curve).U = [];
    for k = Zmod'
        ocvrel.MSMR.(curve).U(end+1, 1) = fun(k);
    end
end
% -------------------------------------------------- Relative lithiation
figure;
for curve = MSMR_curves
    plot(1-ocvrel.estimates(find(string({ocvrel.estimates.name}) == curve)).Z, ...
         ocvrel.estimates(find(string({ocvrel.estimates.name}) == curve)).U); hold on
end
for curve = MSMR_curves % fit
    plot(ocvrel.MSMR.(curve).Zrel, ocvrel.MSMR.(curve).U, 'k:')
end
pptFormatMultiplePlots('Times New Roman')
ylabel('Potential vs. Li (V)');
xlabel('Relative degree of lithiation $\tilde{\theta}_\mathrm{s}$', Interpreter='latex'); 
legend('C/100 dch', 'C/100 chg', 'C/100 avg', 'Model fit', Location='best')
legend(Interpreter="latex", FontName='Times New Roman', FontSize=16, Box='off', Location='best', IconColumnWidth=15)
ylim(v_range); xlim(sol_range)
title('MSMR model fits');

figure;
for curve = MSMR_curves
    idx = find(string({ocvrel.estimates.name}) == curve, 1);
    Zexp = 1 - ocvrel.estimates(idx).Z(:);
    Uexp = ocvrel.estimates(idx).U(:);
    funModel = ocvrel.MSMR.(curve).fun;
    Umodel = zeros(size(Zexp));
    for k = 1:numel(Zexp)
        Umodel(k) = funModel(Zexp(k));
    end
    err = Uexp - Umodel;
    plot(Zexp, err); hold on
end
pptFormatMultiplePlots('Times New Roman');
xlabel('Relative degree of lithiation $\tilde{\theta}_\mathrm{s}$', 'Interpreter', 'latex');
ylabel('$U_\mathrm{MSMR}-U_\mathrm{exp}$ (V)', 'Interpreter', 'latex');
legend('dch', 'chg', 'avg', Location='best')
legend(Interpreter="latex", FontName='Times New Roman', FontSize=16, Box='off', Location='best', IconColumnWidth=15)
ylim([-0.05 0.05]); xlim(sol_range)
title('MSMR fits errors');

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
[Zmod,dZmod,xj,Uj,dZj] = compositionMSMR(ocvrel.MSMR.charge.U, paramsAvg, 25);
Zrel = (Zmod - paramsAvg.thetamin) ./ (paramsAvg.thetamax - paramsAvg.thetamin);
Zrel = min(ocvrel.MSMR.charge.Zrel) + Zrel * (max(ocvrel.MSMR.charge.Zrel) - min(ocvrel.MSMR.charge.Zrel));
U = zeros(size(Zmod));
for k = 1:length(Zmod)
    U(k) = fzero(@(U) lithiationBalanceEquation(U,Zmod(k),paramsAvg,25), mean(ocvrel.MSMR.charge.U));
end
ocvrel.MSMR.average.Zmod = Zmod;
ocvrel.MSMR.average.Zrel = Zrel;
ocvrel.MSMR.average.dZmod = dZmod;
ocvrel.MSMR.average.xj = xj;
ocvrel.MSMR.average.Uj = Uj;
ocvrel.MSMR.average.dZj = dZj;
ocvrel.MSMR.average.U = U;
% -------------------------------------------------- Absolute lithiation
figure
for curve = MSMR_curves
    plot(ocvrel.MSMR.(curve).Zmod, ocvrel.MSMR.(curve).U); hold on
end
pptFormatMultiplePlots('Times New Roman')
ylabel('Potential vs. Li (V)');
xlabel('Absolute degree of lithiation $\tilde{\theta}_\mathrm{s}$', Interpreter='latex');
legend('$U_\mathrm{ocp}^\mathrm{chg}$', '$U_\mathrm{ocp}^\mathrm{dch}$', '$U_\mathrm{ocp}^\mathrm{davg}$', '$U_\mathrm{ocp}^\mathrm{avg}$', Location='best'); 
legend(Interpreter="latex", FontName='Times New Roman', FontSize=16, Box='off', Location='best', IconColumnWidth=15)
ylim(v_range); xlim(sol_range);
title('MSMR OCP models')
drawnow

%% AM correction
Q_MSMR      = ocvrel.Q/(ocvrel.MSMR.average.params.thetamax - ocvrel.MSMR.average.params.thetamin);
AM_per_diff = Q_MSMR/Q_r - AM_per;
AM_per      = AM_per + AM_per_diff;
Q_r_AM      = Q_r*AM_per;                          % corrected sample capacity
ocvrel.theta_max_corr = max(ocvrel.MSMR.average.Zmod);
ocvrel.theta_min_corr = min(ocvrel.MSMR.average.Zmod);
% --------------------------------------------------
idx = find(string({ocvrel.estimates.name}) == "diagonalAverage");
figure
plot(ocvrel.estimates(idx).theta, ocvrel.estimates(idx).U, DisplayName='$U_\mathrm{ocp}^\mathrm{avg, RE}$'); hold on
plot(ocvrel.MSMR.average.Zmod, ocvrel.MSMR.average.U, DisplayName='$U_\mathrm{ocp}^\mathrm{avg, MSMR}$');
pptFormatMultiplePlots('Times New Roman')
xlabel('Absolute degree of lithiation $\theta_\mathrm{s}$', Interpreter='latex');
ylabel('Potential vs. Li (V)');
legend(Interpreter="latex", FontName='Times New Roman', FontSize=16, Box='off', Location='best', IconColumnWidth=15)
ylim(v_range); xlim(sol_range);
title('MSMR vs. Rough estimate of lithiation')
drawnow
% --------------------------------------------------
for k = 1:length(ocvrel.estimates)  % Change absolut lithiation range to match MSMR
    ocvrel.estimates(k).theta = ocvrel.theta_min_corr + (1 - ocvrel.estimates(k).Z) .* (ocvrel.theta_max_corr - ocvrel.theta_min_corr);
end
plot(ocvrel.estimates(idx).theta, ocvrel.estimates(idx).U, DisplayName='$U_\mathrm{ocp}^\mathrm{avg, corr}$');

%% MSMR model plots
figure;
plot(ocvrel.MSMR.average.Zmod, ocvrel.MSMR.average.U, '-',DisplayName='$U_\mathrm{ocp}^\mathrm{MSMR}$'); hold on
pptFormatMultiplePlots('Times New Roman')
plot(ocvrel.MSMR.average.dZmod./max(ocvrel.MSMR.average.dZmod), ocvrel.MSMR.average.U, '-',DisplayName='$d\theta_\mathrm{s}/dU$', LineWidth=1.25);
xlabel('$\theta_\mathrm{s}$ or normalized $d\theta_\mathrm{s}/dU$', Interpreter='latex');
ylabel('Potential vs. Li (V)');
legend(Interpreter="latex", FontName='Times New Roman', FontSize=16, Box='off', Location='best', IconColumnWidth=15)
ylim(v_range); xlim(sol_range);
title('MSMR model')
drawnow
% --------------------------------------------------
figure;
for k = 1:size(ocvrel.MSMR.average.xj, 1)
    plot(ocvrel.MSMR.average.xj(k, :), ocvrel.MSMR.average.Uj(k, :), DisplayName=sprintf('$j$ = %d', k)); hold on
end
pptFormatMultiplePlots('Times New Roman')
xlabel('Fractional occupation of Li $x_\mathrm{j}$', Interpreter='latex');
ylabel('Potential vs. Li (V)');
legend(Interpreter="latex", FontName='Times New Roman', FontSize=16, Box='off', Location='best', IconColumnWidth=15)
title('$U_\mathrm{j}$', Interpreter='latex')
ylim(v_range);
drawnow
% --------------------------------------------------
figure;
for k = 1:size(ocvrel.MSMR.average.dZj, 1)
    plot(ocvrel.MSMR.average.Uj(k, :), ocvrel.MSMR.average.dZj(k, :), DisplayName=sprintf('$j$ = %d', k)); hold on
end
pptFormatMultiplePlots('Times New Roman')
xlabel('Potential vs. Li (V)');
ylabel('$|dx_\mathrm{j}/dU_\mathrm{j}|$', Interpreter='latex');
legend(Interpreter="latex", FontName='Times New Roman', FontSize=16, Box='off', Location='best', IconColumnWidth=15)
xlim(v_range);
drawnow

%% Save model for future use
model.UocpEstimates = ocvrel.estimates;
model.MSMR = ocvrel.MSMR;
% --------------------------------------------------
model.theta_max = ocvrel.theta_max_corr;                                             
model.theta_min = ocvrel.theta_min_corr;
theta_max = model.theta_max;
theta_min = model.theta_min;
model.zToTheta = @(z) z .* (theta_min - theta_max) + theta_max;
model.thetaToZ = @(theta) (theta - theta_max) ./ (theta_min - theta_max);
% --------------------------------------------------
% Z_ch  = model.MSMR.charge.Zmod;
% U_ch  = model.MSMR.charge.U;
% Z_dch = model.MSMR.discharge.Zmod;
% U_dch = model.MSMR.discharge.U;
% Z_avg = model.MSMR.average.Zmod;
% U_avg = model.MSMR.average.U;
Z_ch  = model.UocpEstimates(find(string(OCVestimators) == "charge")).theta;
U_ch  = model.UocpEstimates(find(string(OCVestimators) == "charge")).U;
Z_dch = model.UocpEstimates(find(string(OCVestimators) == "discharge")).theta;
U_dch = model.UocpEstimates(find(string(OCVestimators) == "discharge")).U;
Z_avg = model.UocpEstimates(find(string(OCVestimators) == "diagonalAverage")).theta;
U_avg = model.UocpEstimates(find(string(OCVestimators) == "diagonalAverage")).U;
model.UocpCh = @(theta) interp1(Z_ch, U_ch, theta, 'linear', 'extrap');
model.UocpDch = @(theta) interp1(Z_dch, U_dch, theta, 'linear', 'extrap');
model.UocpAvg = @(theta) interp1(Z_avg, U_avg, theta, 'linear', 'extrap');
% --------------------------------------------------
model.Q25 = ocvrel.Q;
model.Qnom = Q_r_AM;
% --------------------------------------------------
fileName = sprintf('output/Uocp_n.mat');
save(fileName, 'model');