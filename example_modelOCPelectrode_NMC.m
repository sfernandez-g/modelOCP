%EXAMPLE_MODELOCPELECTRODE_NMC Fit an NMC811 half-cell OCP model at 25 degC.
% Run this script from the MATLAB editor or with:
%   run('example_modelOCPelectrode_NMC.m')
% Requires the experimental workbook, capacity MAT file, and fitting toolboxes.

repoRoot = fileparts(mfilename('fullpath'));
addpath(fullfile(repoRoot, 'fun'));
settings = struct;

%% Experimental inputs
settings.inputFile = fullfile(repoRoot, 'input', 'pos', ...
    'checkUp_LGMH1_PAT20261302_NMC3_Li_T25_Channel_16_Wb_1.xlsx');
settings.capacityFile = fullfile(repoRoot, 'output', 'theoPracCapPos.mat'); % Q_r [mAh].
settings.cellName = 'LG MH1 PAT half-cell PE';
settings.electrodeName = 'NMC811';
settings.temperatureDegC = 25;
settings.voltageRange = [3 4.3]; % Potential versus Li [V].
settings.activeMaterialFraction = 0.96; % Assumed NMC active-material fraction.
settings.lithiationAnchor = "lithiated"; % Assume theta_max = 1.

%% MSMR fitting inputs
settings.initialMSMR = []; % Estimate initial parameters from the measured data.
settings.galleryCount = 4;
settings.boundScales = [0.3 0.2 0.1]; % Lower-bound factors for U0, X, w.
settings.differentialWeight = 1e-2;
settings.optimizer.seed = 0; % Reproducible GA initialization.
settings.hiddenGalleryCorrection = true; % Apply the assumed hidden-gallery correction.

%% Optional figure settings
% Each inset row is [left bottom width height xmin xmax].
% Positions are normalized; x limits are dimensionless lithiation.
settings.relativeInsets = [.7 .35 .125 .125 .975 1; ...
                          .2 .6 .125 .125 0 .01];
settings.fitInsets = [.6 .75 .125 .125 .95 1];
settings.absoluteInsets = [.6 .75 .125 .125 .95 1];
settings.showRelativeLegend = false;
settings.showComparisonLegend = false;

%% Run the shared workflow
settings.outputFile = ""; % Return results without saving a model MAT file.
% To save the model, set an output path before the function call, for example:
% settings.outputFile = fullfile(repoRoot, 'output', 'Uocp_p_example.mat');
[model, ocvrel, supplementaryResults, ocpFigures] = modelOCPelectrode(settings);

% Returned results:
% model                Final electrode model, including MSMR and LUT fields.
% ocvrel               Processed measurements and fitted curves.
% supplementaryResults Snapshot before hidden-gallery and AM corrections.
% ocpFigures           Seven figure handles, indexed in window-number order.
