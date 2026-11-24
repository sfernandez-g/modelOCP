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
file      = fullfile(repoRoot, 'input', 'commercial cell', 'checkUp_paramMW_LGMH1_T25_Channel_3_Wb_1.xlsx');
%--------------------------------------------------------------------------
cellName  = 'LG MH1 commercial cell';
%--------------------------------------------------------------------------
v_range   = [2.5 4.2];
% -------------------------------------------------------------------------
sol_range   = [0 1];
% -------------------------------------------------------------------------
test_data = getDataArbin_XLSX(file, 7); % filename, step where TC increments
% -------------------------------------------------------------------------
figure('Name', [cellName, ': measured voltage'], 'NumberTitle', 'off');
plot(test_data.t./3600, test_data.v, DisplayName='Voltage'); hold on;
xlabel('Time (h)'); ylabel('Voltage (V)'); legend(Location="best");
thesisFormat;
% -------------------------------------------------------------------------
cfg = struct('name', {'CC_dis','CV_dis','CC_chg','CV_chg'}, ... % test phases
             'step', {5, 6, 3, 4}, ...                          % step values of phases for processing
             'counter', {1, 1, 1, 1} ...                        % counter values of phases for processing
             );
voltChannels = {'MainVolt'};
filtv        = [2.5, 4.2];
[cal, ocvrel, OCVestimators] = processOCPCurves(test_data, cfg, voltChannels, 1, cellName, filtv, @(Z) 0 < Z & Z < 0.1, 1e4);
% -------------------------------------------------------------------------
figure('Name', [cellName, ': OCV versus SOC'], 'NumberTitle', 'off');
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
model.UocpEstimates = ocvrel.estimates;
% ---------------------------------------------
Z_ch  = model.UocpEstimates(find(string(OCVestimators) == "charge")).Z;
U_ch  = model.UocpEstimates(find(string(OCVestimators) == "charge")).U;
Z_dch = model.UocpEstimates(find(string(OCVestimators) == "discharge")).Z;
U_dch = model.UocpEstimates(find(string(OCVestimators) == "discharge")).U;
Z_avg = model.UocpEstimates(find(string(OCVestimators) == "diagonalAverage")).Z;
U_avg = model.UocpEstimates(find(string(OCVestimators) == "diagonalAverage")).U;
model.UocvCh = @(theta) interp1(Z_ch, U_ch, theta, 'linear', 'extrap');
model.UocvDch = @(theta) interp1(Z_dch, U_dch, theta, 'linear', 'extrap');
model.UocvAvg = @(theta) interp1(Z_avg, U_avg, theta, 'linear', 'extrap');
model.Q25 = ocvrel.Q;

fileName = fullfile(repoRoot, 'output', 'Uocv_comCell.mat');
if saveModel
    save(fileName, 'model');
end
