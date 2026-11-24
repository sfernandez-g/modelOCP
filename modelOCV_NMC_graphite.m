clc; clear; close all;
repoRoot = fileparts(mfilename('fullpath'));
RGB = orderedcolors("gem");
saveModel = false; % Set true to explicitly export the processed model.
if saveModel && ~isfolder(fullfile(repoRoot, 'output'))
    mkdir(fullfile(repoRoot, 'output'));
end

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
figure('Name', [cellName, ': measured voltage'], 'NumberTitle', 'off');
plot(test_data.t./3600, test_data.v, DisplayName='Voltage'); hold on;
plot(test_data.t./3600, test_data.Aux_V1, DisplayName='AuxV');
plot(test_data.t./3600, test_data.Aux_V1-test_data.v, DisplayName='AuxV - Voltage');
xlabel('Time (h)'); ylabel('Voltage (V)'); legend(Location="best");
thesisFormat;
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
    figure('Name', sprintf('%s: %s versus SOC', cellName, voltChannels{j}), 'NumberTitle', 'off');
    for k = 1:length(ocvrel.estimates)
        plot(ocvrel.estimates(k).Z, ocvrel.estimates(k).U); hold on
    end
    xlabel('SOC $z$', Interpreter='latex');
    ylabel('Voltage or Potential vs. Li (V)');
    legend('$U^\mathrm{dch}$', '$U^\mathrm{chg}$', '$U^\mathrm{avg}$', Location='best');
    xlim(sol_range);
    thesisFormat;
    drawnow

    %% Assemble the model; optional export
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
    if saveModel
        save(fileName, 'model');
    end
end
