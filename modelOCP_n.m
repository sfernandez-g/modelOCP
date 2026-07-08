clc; clear; close all; RGB = orderedcolors("gem");
addpath("fun\"); addpath("fun\Uoc\");

%% Sample info
file        = "input\neg\checkUp_LGMH1_PAT20261302_G1_Li_T25_2026_02_20_155142\checkUp_LGMH1_PAT20261302_G1_Li_T25_Channel_2_Wb_1.xlsx";
cellName    = 'LG MH1 half-cell neg';
v_range     = [0 1.5];
sol_range   = [0 1];

%% Sample capacity
load("output\theoPracCapNeg.mat");
Q_r     = Q_r*1e-3;
AM_per  = 0.98;          % Rojas 2024
Q_r_AM  = Q_r*AM_per

%% Test data
test_data   = getDataArbin_XLSX(file, 7); % filename, step where TC increments
figure; plot(test_data.t./3600, test_data.v);
articleFormatMultiplePlots('Times New Roman')
xlabel('Time (h)', Interpreter='latex'); 
ylabel('Potential vs. Li (V)', Interpreter='latex'); 
ylim(v_range);
title('Check-up test');

%% OCP-Theta_rel estimation
cfg                 = struct('name', {'CC_dis','CV_dis','CC_chg','CV_chg'}, ... % test phases
                             'step', {5, 6, 3, 4}, ...                          % step values of phases for processing
                             'counter', {1, 1, 1, 1} ...                        % counter values of phases for processing
                             );
cal                 = processCU(test_data, cfg);
OCVresolution       = 1e4; % this many points between 0 and 1
OCVestimators       = {'diagonalAverage'};
cellData.raw        = cal; 
cellData.cal        = cal; 
cellData.name       = cellName;
config.vmin         = cal.minV; config.vmax = cal.maxV;
config.datype       = 'useAvg';
config.daxcorrfiltz = @(Z){0<Z & Z<0.1};   % focus on this SOC region
config.daxcorrfiltv = @(V){0<V & V<1.5};   % focus on this V region
ocvrel              = estimateOCV(cellData, config, OCVresolution, OCVestimators);
ocvrel                  = addDchChOCV(ocvrel, cellData, "discharge");
ocvrel                  = addDchChOCV(ocvrel, cellData, "charge");
OCVestimators{end+1}    = 'charge';
OCVestimators{end+1}    = 'discharge';
ocvrel.estimates        = ocvrel.estimates(end:-1:1);
ocvrel.Q = cal.QAh

figure;
for k = 1:length(ocvrel.estimates)
    plot(1-ocvrel.estimates(k).Z, ocvrel.estimates(k).U); hold on
end
articleFormatMultiplePlots('Times New Roman')
xlabel('Relative degree of lithiation $\tilde{\theta}_\mathrm{s}$', Interpreter='latex'); 
ylabel('Potential vs. Li (V)');
legend('C/30 chg', 'C/30 dch', 'C/30 avg', Location='best'); 
legend(Interpreter="latex", FontName='Times New Roman', FontSize=8, Box='off', Location='best', IconColumnWidth=15)
ylim(v_range); xlim(sol_range);
title('OCP vs. Relative degree of lithiation');
drawnow

%% Rough estimate of absolute lithiation
ocvrel.theta_max = max((1-ocvrel.estimates(2).Z)).*ocvrel.Q/Q_r_AM; % Based on "discharge" data asumming it starts from a completely delithiated state
ocvrel.theta_min = 0;                                               % Assume we reach minimum lithiation

for k = 1:length(ocvrel.estimates)                                         
    Z = 1 - ocvrel.estimates(k).Z;                                  % Change SOC to relative lithiation
    ocvrel.estimates(k).theta = ocvrel.theta_min + Z .* (ocvrel.theta_max - ocvrel.theta_min);
end

figure;
for k = 1:length(ocvrel.estimates)
    plot(ocvrel.estimates(k).theta, ocvrel.estimates(k).U)
    hold on
end
articleFormatMultiplePlots('Times New Roman')
xlabel('Absolute degree of lithiation $\theta_\mathrm{s}$', Interpreter='latex'); 
ylabel('Potential vs. Li (V)', Interpreter='latex');
legend('C/30 chg', 'C/30 dch', 'C/30 avg', Location='best'); 
legend(Interpreter="latex", FontName='Times New Roman', FontSize=8, Box='off', Location='best', IconColumnWidth=15)
ylim(v_range); xlim(sol_range);
title('OCP vs. Absolute degree of lithiation (rough est.)');
drawnow

%% MSMR model fit
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
for curve = MSMR_curves % data
    plot(1-ocvrel.estimates(find(string({ocvrel.estimates.name}) == curve)).Z, ...
         ocvrel.estimates(find(string({ocvrel.estimates.name}) == curve)).U); hold on
end
for curve = MSMR_curves % fit
    plot(ocvrel.MSMR.(curve).Zrel, ocvrel.MSMR.(curve).U, 'k:')
end
articleFormatMultiplePlots('Times New Roman')
ylabel('Potential vs. Li (V)', Interpreter='latex');
xlabel('Relative degree of lithiation $\tilde{\theta}_\mathrm{s}$', Interpreter='latex'); 
legend('C/100 dch', 'C/100 chg', 'C/100 avg', 'Model fit', '', '',Location='best')
legend(Interpreter="latex", FontName='Times New Roman', FontSize=8, Box='off', Location='best', IconColumnWidth=15)
ylim(v_range); xlim(sol_range)
title('MSMR model fits');

% -------------------------------------------------- Average MSMR
% fields = fieldnames(params);
% for i = 1:length(fields)
%     ocvrel.MSMR.average.params.(fields{i}) = mean([ocvrel.MSMR.discharge.params.(fields{i}), ocvrel.MSMR.charge.params.(fields{i})], 2);
% end
% ocvrel.MSMR.average.Zmod = ocvrel.MSMR.diagonalAverage.Zmod;
% ocvrel.MSMR.average.fun = @(z) fzero(@(U) lithiationBalanceEquation(U, z, ocvrel.MSMR.average.params, 25), mean(ocvrel.MSMR.diagonalAverage.U));
% ocvrel.MSMR.average.U = [];
% for k = ocvrel.MSMR.average.Zmod'
%     ocvrel.MSMR.average.U(end+1, 1) = ocvrel.MSMR.average.fun(k);
% end
% 
% [ocvrel.MSMR.average.Zmod, ocvrel.MSMR.average.dZmod, ...
%  ocvrel.MSMR.average.xj, ocvrel.MSMR.average.Uj, ... 
%  ocvrel.MSMR.average.dZj] = compositionMSMR(ocvrel.MSMR.average.U, ocvrel.MSMR.average.params, 25);
% 
% U_new = linspace(v_range(1), v_range(2), 10000);
% ocvrel.MSMR.average.Zmod = interp1(ocvrel.MSMR.average.U, ocvrel.MSMR.average.Zmod, U_new, "linear", "extrap");
% ocvrel.MSMR.average.Zrel = (ocvrel.MSMR.average.Zmod - ocvrel.MSMR.average.params.thetamin)./...
%                            (ocvrel.MSMR.average.params.thetamax - ocvrel.MSMR.average.params.thetamin);
% ocvrel.MSMR.average.Zrel = min(ocvrel.MSMR.diagonalAverage.Zrel) + ocvrel.MSMR.average.Zrel.*(max(ocvrel.MSMR.diagonalAverage.Zrel) - min(ocvrel.MSMR.diagonalAverage.Zrel));
% ocvrel.MSMR.average.dZmod = interp1(ocvrel.MSMR.average.U, ocvrel.MSMR.average.dZmod, U_new, "linear", "extrap");
% ocvrel.MSMR.average.U = U_new;
% 
% plot(ax, ocvrel.MSMR.average.Zrel, ocvrel.MSMR.average.U, LineWidth=1.5)

% -------------------------------------------------- Inset ax
inset_ax = axes('Position', [.25, .6, .2, .25], 'Box', 'on', 'FontSize', 8); 
for curve = MSMR_curves
    plot(inset_ax, ...
         1 - ocvrel.estimates(find(string({ocvrel.estimates.name}) == curve)).Z, ...
         ocvrel.estimates(find(string({ocvrel.estimates.name}) == curve)).U, LineWidth=1.5); 
    hold(inset_ax, 'on');
end
for curve = MSMR_curves
    plot(inset_ax, ocvrel.MSMR.(curve).Zrel, ocvrel.MSMR.(curve).U, 'k:', LineWidth=1.5);
end
% plot(inset_ax, ocvrel.MSMR.average.Zrel, ocvrel.MSMR.average.U, LineWidth=1.5)
xlim(inset_ax, [0.05, 0.15]);
set(inset_ax, 'FontName', 'Times New Roman');

inset_ax = axes('Position', [.625, .5, .2, .25], 'Box', 'on', 'FontSize', 8); 
for curve = MSMR_curves
    plot(inset_ax, ...
         1 - ocvrel.estimates(find(string({ocvrel.estimates.name}) == curve)).Z, ...
         ocvrel.estimates(find(string({ocvrel.estimates.name}) == curve)).U, LineWidth=1.5); 
    hold(inset_ax, 'on');
end
for curve = MSMR_curves
    plot(inset_ax, ocvrel.MSMR.(curve).Zrel, ocvrel.MSMR.(curve).U, 'k:', LineWidth=1.5);
end
% plot(inset_ax, ocvrel.MSMR.average.Zrel, ocvrel.MSMR.average.U, LineWidth=1.5)
xlim(inset_ax, [0.45, 0.6]);
set(inset_ax, 'FontName', 'Times New Roman');
drawnow

% -------------------------------------------------- Absolute lithiation
MSMR_curves = ["discharge", "charge", "diagonalAverage"];
figure
for curve = MSMR_curves % fit
    plot(ocvrel.MSMR.(curve).Zmod, ocvrel.MSMR.(curve).U); hold on
end
articleFormatMultiplePlots('Times New Roman')
ylabel('Potential vs. Li (V)');
xlabel('Absolute degree of lithiation $\tilde{\theta}_\mathrm{s}$', Interpreter='latex');
legend('$U_\mathrm{ocp}^\mathrm{chg}$', '$U_\mathrm{ocp}^\mathrm{dch}$', '$U_\mathrm{ocp}^\mathrm{avg}$', Location='best'); 
legend(Interpreter="latex", FontName='Times New Roman', FontSize=8, Box='off', Location='best', IconColumnWidth=15)
ylim(v_range); xlim(sol_range);
title('MSMR OCP models')

% -------------------------------------------------- AM correction
Q_MSMR      = ocvrel.Q/(ocvrel.MSMR.diagonalAverage.params.thetamax - ocvrel.MSMR.diagonalAverage.params.thetamin);
AM_per_diff = Q_MSMR/Q_r - AM_per;
AM_per      = AM_per + AM_per_diff
Q_r_AM      = Q_r*AM_per                            % corrected sample capacity

ocvrel.theta_max_corr = max(ocvrel.MSMR.diagonalAverage.Zmod);
ocvrel.theta_min_corr = min(ocvrel.MSMR.diagonalAverage.Zmod);

idx = find(string({ocvrel.estimates.name}) == "diagonalAverage");
figure
plot(ocvrel.estimates(idx).theta, ocvrel.estimates(idx).U, DisplayName='$U_\mathrm{ocp}^\mathrm{avg, RE}$'); hold on
plot(ocvrel.MSMR.diagonalAverage.Zmod, ocvrel.MSMR.diagonalAverage.U, DisplayName='$U_\mathrm{ocp}^\mathrm{avg, MSMR}$');
articleFormatMultiplePlots('Times New Roman')
xlabel('Absolute degree of lithiation $\theta_\mathrm{s}$', Interpreter='latex');
ylabel('Potential vs. Li (V)');
legend(Interpreter="latex", FontName='Times New Roman', FontSize=8, Box='off', Location='best', IconColumnWidth=15)
ylim(v_range); xlim(sol_range);
title('MSMR vs. Rough estimate of lithiation')

% -------------------------------------------------- Final MSMR model
figure;
plot(ocvrel.MSMR.diagonalAverage.Zmod, ocvrel.MSMR.diagonalAverage.U, '-',DisplayName='$U_\mathrm{ocp}^\mathrm{MSMR}$'); hold on
articleFormatMultiplePlots('Times New Roman')
plot(ocvrel.MSMR.diagonalAverage.dZmod./max(ocvrel.MSMR.diagonalAverage.dZmod), ocvrel.MSMR.diagonalAverage.U, '-',DisplayName='$d\theta_\mathrm{s}/dU$', LineWidth=1.25);
xlabel('$\theta_\mathrm{s}$ or normalized $d\theta_\mathrm{s}/dU$', Interpreter='latex');
ylabel('Potential vs. Li (V)');
legend(Interpreter="latex", FontName='Times New Roman', FontSize=8, Box='off', Location='best', IconColumnWidth=15)
ylim(v_range); xlim(sol_range);
title('MSMR model')

figure;
for k = 1:size(ocvrel.MSMR.diagonalAverage.xj, 1)
    plot(ocvrel.MSMR.diagonalAverage.xj(k, :), ocvrel.MSMR.diagonalAverage.Uj(k, :), DisplayName=sprintf('$j$ = %d', k)); hold on
end
articleFormatMultiplePlots('Times New Roman')
xlabel('Fractional occupation of Li $x_\mathrm{j}$', Interpreter='latex');
ylabel('Potential vs. Li (V)');
legend(Interpreter="latex", FontName='Times New Roman', FontSize=8, Box='off', Location='best', IconColumnWidth=15)
title('$U_\mathrm{j}$', Interpreter='latex')
ylim(v_range);

figure;
for k = 1:size(ocvrel.MSMR.diagonalAverage.dZj, 1)
    plot(ocvrel.MSMR.diagonalAverage.Uj(k, :), ocvrel.MSMR.diagonalAverage.dZj(k, :), DisplayName=sprintf('$j$ = %d', k)); hold on
end
articleFormatMultiplePlots('Times New Roman')
xlabel('Potential vs. Li (V)');
ylabel('$|dx_\mathrm{j}/dU_\mathrm{j}|$', Interpreter='latex');
legend(Interpreter="latex", FontName='Times New Roman', FontSize=8, Box='off', Location='best', IconColumnWidth=15)
xlim(v_range);

%% Save model for future use
model.UocpEstimates = ocvrel.estimates;
model.Q25 = ocvrel.Q;
model.Qnom = Q_r_AM;
model.theta_max = ocvrel.theta_max_corr;                                             
model.theta_min = ocvrel.theta_min_corr;
model.MSMR = ocvrel.MSMR;

model.zToTheta = @(z) z .* (model.theta_min - model.theta_max) + model.theta_max;
model.thetaToZ = @(theta) (theta - model.theta_max) ./ (model.theta_min - model.theta_max);
model.UocpCh = @(theta) interp1(model.MSMR.charge.Zmod, ...
                                model.MSMR.charge.U, ...
                                theta, 'linear', 'extrap');
model.UocpDch = @(theta) interp1(model.MSMR.discharge.Zmod, ...
                                 model.MSMR.discharge.U, ...
                                 theta, 'linear', 'extrap');
model.UocpAvg = @(theta) interp1(model.MSMR.diagonalAverage.Zmod, ...
                                 model.MSMR.diagonalAverage.U, ...
                                 theta, 'linear', 'extrap');

fileName = sprintf('output/Uocp_n.mat');
save(fileName, 'model');