clc; clear; close all; 
RGB = orderedcolors("gem");

%% add utility paths
addpath("fun\"); addpath("fun\Uoc\");
%% process OCV data
file      = 'input\full-cell\CU_qOCP_LGMH1_PAT20251015_NMC1_G5_Channel_9_Wb_1.xlsx';
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
    cal = processCU(test_data, cfg, voltChannels, j); % test data, processing config, voltChannels
    %% Estimate OCP
    cellData.raw = cal; cellData.cal = cal; cellData.name = cellName;
    % -------------------------------------------------------------------------
    config.vmin         = cal.minV; config.vmax = cal.maxV; config.datype = 'useAvg';
    config.daxcorrfiltz = @(Z){0<Z & Z<0.1};                                            % focus on this SOC region
    config.daxcorrfiltv = @(V){filtv(j, 1)<V & V<filtv(j, 2)};                          % focus on this V region
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

    fileName = sprintf('output/Uocv_%s.mat', voltChannels{j});
    save(fileName, 'model');
end