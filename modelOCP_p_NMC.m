clc; clear; close all; 
RGB = orderedcolors("gem");

%% add utility paths
addpath("fun\"); addpath("fun\Uoc\");

%% process OCP data
file      = 'data\pos\checkUp_LGMH1_PAT20261302_NMC3_Li_T25_2026_03_02_094138\checkUp_LGMH1_PAT20261302_NMC3_Li_T25_Channel_16_Wb_1.xlsx';
%--------------------------------------------------------------------------
cellName  = 'LG MH1 PAT half-cell PE'; 
%--------------------------------------------------------------------------
v_range   = [3 4.3];
% -------------------------------------------------------------------------
sol_range   = [0 1];
% -------------------------------------------------------------------------
load("output\theoPracCapPos.mat");
Q_r     = Q_r*1e-3;
AM_per  = 0.96;          % Rojas 2024
Q_r_AM  = Q_r*AM_per;
% -------------------------------------------------------------------------
test_data = getDataArbin_XLSX(file); % filename, step where TC increments
% -------------------------------------------------------------------------
figure;
plot(test_data.t./3600, test_data.v, DisplayName='Voltage'); hold on;
pptFormatMultiplePlots('Times New Roman')
xlabel('Time (h)'); ylabel('Voltage (V)'); legend(Location="best");
% -------------------------------------------------------------------------
cfg = struct('name', {'CC_dis','CV_dis','CC_chg','CV_chg'}, ... % test phases
             'step', {5, 6, 3, 4}, ...                          % step values of phases for processing
             'counter', {0, 0, 1, 1} ...                        % counter values of phases for processing
             );
voltChannels = {'MainVolt'};
filtv        = [3, 4.3];
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

%% Rough estimate of absolute lithiation
ocvrel.theta_max = 1;                                                      % Assume we reach maximum lithiation
ocvrel.theta_min = 1 - max((1-ocvrel.estimates(1).Z)).*ocvrel.Q/Q_r_AM;    % Based on "charge" data asuming it starts from a completely lithiated state                                             % Assume we reach minimum lithiation
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
config.init.U0 = [3.580811; 3.739531; 3.998416; 4.18913];  % Garrick 2023
config.init.X = [0.164246; 0.367062; 0.31328; 0.155412];
config.init.w = [1.281307; 1.918194; 2.714622; 0.787826];
config.init.thetamin = 0;
config.init.thetamax = 1;
% ---------------------------------------------------
len = length(config.init.U0);
% -------------------------------------------------- Literature based bounds
config.lb.U0 = config.init.U0.*0.5;
config.lb.X = config.init.X.*0.5;
config.lb.w = config.init.w.*0.3;
config.lb.thetamin = 0-1e-10;
config.lb.thetamax = 0.85;

config.ub.U0 = config.init.U0./0.5;
config.ub.X = config.init.X./0.5;
config.ub.w = config.init.w./0.3;
config.ub.thetamin = 0.1;
config.ub.thetamax = 1+1e-10;
% -------------------------------------------------- Aditional bounds or constraints
A = []; % n x nvars
b = [];  % n x 1
cons.A = A;
cons.b = b;
% -------------------------------------------------- Differential curve weight
config.optw = 1e-10;
% -------------------------------------------------- Launch MSMR fit
vis = 0;
MSMR_curves = ["diagonalAverage", "charge", "discharge"];
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
    config.init.U0 = ocvrel.MSMR.diagonalAverage.params.U0;
    config.init.X = ocvrel.MSMR.diagonalAverage.params.X;
    config.init.w = ocvrel.MSMR.diagonalAverage.params.w;
    config.init.thetamin = ocvrel.MSMR.diagonalAverage.params.thetamin;
    config.init.thetamax = ocvrel.MSMR.diagonalAverage.params.thetamax;

    config.lb.U0 = config.init.U0.*0.8;
    config.lb.X = config.init.X.*0.6;
    config.lb.w = config.init.w.*0.4;
    config.lb.thetamin = config.init.thetamin - 1e-2;
    config.lb.thetamax = config.init.thetamax - 1e-2;
    
    config.ub.U0 = config.init.U0./0.8;
    config.ub.X = config.init.X./0.6;
    config.ub.w = config.init.w./0.4;
    config.ub.thetamin = config.init.thetamin + 1e-2;
    config.ub.thetamax = config.init.thetamax + 1e-2;
end
% -------------------------------------------------- Plot Fits
MSMR_curves = ["discharge", "charge", "diagonalAverage"];
figure; ax = gca;
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
legend('C/100 dch', 'C/100 chg', 'Model fit', '', Location='best')
legend(Interpreter="latex", FontName='Times New Roman', FontSize=16, Box='off', Location='best', IconColumnWidth=15)
ylim(v_range); xlim(sol_range)
title('MSMR model fits');
% -------------------------------------------------- Absolute lithiation
figure
for curve = MSMR_curves
    plot(ocvrel.MSMR.(curve).Zmod, ocvrel.MSMR.(curve).U); hold on
end
pptFormatMultiplePlots('Times New Roman')
ylabel('Potential vs. Li (V)');
xlabel('Absolute degree of lithiation $\tilde{\theta}_\mathrm{s}$', Interpreter='latex');
legend('$U_\mathrm{ocp}^\mathrm{chg}$', '$U_\mathrm{ocp}^\mathrm{dch}$', '$U_\mathrm{ocp}^\mathrm{avg}$', Location='best'); 
legend(Interpreter="latex", FontName='Times New Roman', FontSize=16, Box='off', Location='best', IconColumnWidth=15)
ylim(v_range); xlim(sol_range);
title('MSMR OCP models')
drawnow

%% Hidden gallery problem
for curve = MSMR_curves
    ocvrel.MSMR.(curve).params.thetamin = ocvrel.theta_min;
    ocvrel.MSMR.(curve).params.thetamax = ocvrel.MSMR.(curve).params.thetamax;
    ocvrel.MSMR.(curve).Zmod = ocvrel.MSMR.(curve).Zmod .* (ocvrel.MSMR.(curve).params.thetamax - ocvrel.MSMR.(curve).params.thetamin) + ...
                               ocvrel.MSMR.(curve).params.thetamin;
end
ocvrel.theta_max_corr = max(ocvrel.MSMR.diagonalAverage.Zmod);
ocvrel.theta_min_corr = min(ocvrel.MSMR.diagonalAverage.Zmod);
% -------------------------------------------------- AM correction
Q_MSMR      = ocvrel.Q/(ocvrel.MSMR.diagonalAverage.params.thetamax - ocvrel.MSMR.diagonalAverage.params.thetamin);
AM_per_diff = Q_MSMR/Q_r - AM_per;
AM_per      = AM_per + AM_per_diff;
Q_r_AM      = Q_r*AM_per;                             % corrected sample capacity
% --------------------------------------------------
idx = find(string({ocvrel.estimates.name}) == "diagonalAverage");
figure
plot(ocvrel.estimates(idx).theta, ocvrel.estimates(idx).U, DisplayName='$U_\mathrm{ocp}^\mathrm{avg, RE}$'); hold on
plot(ocvrel.MSMR.diagonalAverage.Zmod, ocvrel.MSMR.diagonalAverage.U, DisplayName='$U_\mathrm{ocp}^\mathrm{avg, MSMR}$');
pptFormatMultiplePlots('Times New Roman')
xlabel('Absolute degree of lithiation $\theta_\mathrm{s}$', Interpreter='latex');
ylabel('Potential vs. Li (V)');
legend(Interpreter="latex", FontName='Times New Roman', FontSize=16, Box='off', Location='best', IconColumnWidth=15)
ylim(v_range); xlim(sol_range);
title('MSMR vs. Rough estimate of lithiation')
drawnow

%% MSMR model plots
figure;
plot(ocvrel.MSMR.diagonalAverage.Zmod, ocvrel.MSMR.diagonalAverage.U, '-',DisplayName='$U_\mathrm{ocp}^\mathrm{MSMR}$'); hold on
pptFormatMultiplePlots('Times New Roman')
plot(ocvrel.MSMR.diagonalAverage.dZmod./max(ocvrel.MSMR.diagonalAverage.dZmod), ocvrel.MSMR.diagonalAverage.U, '-',DisplayName='$d\theta_\mathrm{s}/dU$', LineWidth=1.25);
xlabel('$\theta_\mathrm{s}$ or normalized $d\theta_\mathrm{s}/dU$', Interpreter='latex');
ylabel('Potential vs. Li (V)');
legend(Interpreter="latex", FontName='Times New Roman', FontSize=16, Box='off', Location='best', IconColumnWidth=15)
ylim(v_range); xlim(sol_range);
title('MSMR model')
drawnow
% --------------------------------------------------
figure;
for k = 1:size(ocvrel.MSMR.diagonalAverage.xj, 1)
    plot(ocvrel.MSMR.diagonalAverage.xj(k, :), ocvrel.MSMR.diagonalAverage.Uj(k, :), DisplayName=sprintf('$j$ = %d', k)); hold on
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
for k = 1:size(ocvrel.MSMR.diagonalAverage.dZj, 1)
    plot(ocvrel.MSMR.diagonalAverage.Uj(k, :), ocvrel.MSMR.diagonalAverage.dZj(k, :), DisplayName=sprintf('$j$ = %d', k)); hold on
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
Z_ch  = model.MSMR.charge.Zmod;
U_ch  = model.MSMR.charge.U;
Z_dch = model.MSMR.discharge.Zmod;
U_dch = model.MSMR.discharge.U;
Z_avg = model.MSMR.diagonalAverage.Zmod;
U_avg = model.MSMR.diagonalAverage.U;
model.UocpCh = @(theta) interp1(Z_ch, U_ch, theta, 'linear', 'extrap');
model.UocpDch = @(theta) interp1(Z_dch, U_dch, theta, 'linear', 'extrap');
model.UocpAvg = @(theta) interp1(Z_avg, U_avg, theta, 'linear', 'extrap');
% --------------------------------------------------
model.Q25 = ocvrel.Q;
model.Qnom = Q_r_AM;
% --------------------------------------------------
fileName = sprintf('results/Uocp_p.mat');
save(fileName, 'model');
