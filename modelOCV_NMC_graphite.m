clc; clear; close all; 
repoRoot = fileparts(mfilename('fullpath'));
RGB = orderedcolors("gem");

%% add utility paths
addpath(fullfile(repoRoot, "fun"));

%% process OCV data
file      = fullfile(repoRoot, 'input', 'full-cell', 'CU_qOCP_LGMH1_PAT20251015_NMC1_G5_Channel_9_Wb_1.xlsx');
%--------------------------------------------------------------------------
cellName  = 'LG MH1 PAT full-cell'; 
%--------------------------------------------------------------------------
v_range   = [2.5 4.2];
% -------------------------------------------------------------------------
sol_range   = [0 1];
% -------------------------------------------------------------------------
test_data = getDataArbin_XLSX(file, 7); % filename, step where TC increments
% -------------------------------------------------------------------------
figure;
plot(test_data.t./3600, test_data.v, DisplayName='Voltage'); hold on;
plot(test_data.t./3600, test_data.Aux_V1, DisplayName='AuxV');
plot(test_data.t./3600, test_data.Aux_V1-test_data.v, DisplayName='AuxV - Voltage');
pptFormatMultiplePlots('Times New Roman')
xlabel('Time (h)'); ylabel('Voltage (V)'); legend(Location="best");
% -------------------------------------------------------------------------
cfg = struct('name', {'CC_dis','CV_dis','CC_chg','CV_chg'}, ... % test phases
             'step', {5, 6, 3, 4}, ...                          % step values of phases for processing
             'counter', {1, 1, 1, 1} ...                        % counter values of phases for processing
             );
voltChannels = {'MainVolt', 'AuxVolt', 'AuxMinusMain'};
filtv        = [2.5, 4.2; 3, 4.3; 0, 1.5];

for j = 1:length(voltChannels)
    [cal, ocvrel, OCVestimators] = processOCPCurves(test_data, cfg, voltChannels, j, cellName, filtv(j,:), @(Z) 0 < Z & Z < 0.1, 1e4);
    % -------------------------------------------------------------------------
    figure;
    for k = 1:length(ocvrel.estimates)
        plot(ocvrel.estimates(k).Z, ocvrel.estimates(k).U); hold on
    end
    pptFormatMultiplePlots('Times New Roman')
    xlabel('SOC $z$', Interpreter='latex'); 
    ylabel('Voltage or Potential vs. Li (V)');
    legend('$U^\mathrm{chg}$', '$U^\mathrm{dch}$', '$U^\mathrm{avg}$', Location='best'); 
    legend(Interpreter="latex", FontName='Times New Roman', FontSize=16, Box='off', Location='best', IconColumnWidth=15)
    xlim(sol_range);
    title('OCP vs. Relative degree of lithiation');
    drawnow

    %% Save model for future use
    % ---- Uocp estimates
    model.raw = ocvrel.estimates;
    % ---- LUT
    model.LUT.Z_ch  = model.raw(find(string(OCVestimators) == "charge")).Z;
    model.LUT.U_ch  = model.raw(find(string(OCVestimators) == "charge")).U;
    model.LUT.UocpCh = @(z) linearinterp(model.LUT.Z_ch, model.LUT.U_ch, z);

    model.LUT.Z_dch = model.raw(find(string(OCVestimators) == "discharge")).Z;
    model.LUT.U_dch = model.raw(find(string(OCVestimators) == "discharge")).U;
    model.LUT.UocpDch = @(z) linearinterp(model.LUT.Z_dch, model.LUT.U_dch, z);

    model.LUT.Z_avg = model.raw(find(string(OCVestimators) == "diagonalAverage")).Z;
    model.LUT.U_avg = model.raw(find(string(OCVestimators) == "diagonalAverage")).U;
    model.LUT.UocpAvg = @(z) linearinterp(model.LUT.Z_avg, model.LUT.U_avg, z);
    % --------------------------------------------------
    model.Q25 = ocvrel.Q;

    fileName = fullfile(repoRoot, 'output', sprintf('Uocv_%s.mat', voltChannels{j}));
    save(fileName, 'model');
end
