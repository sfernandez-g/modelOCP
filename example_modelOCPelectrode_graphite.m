%EXAMPLE_MODELOCPELECTRODE_GRAPHITE Fit a graphite half-cell OCP model at 25 degC.
% Run this script from the MATLAB editor or with:
%   run('example_modelOCPelectrode_graphite.m')
% Requires the experimental workbook, capacity MAT file, and fitting toolboxes.

repoRoot = fileparts(mfilename('fullpath'));
addpath(fullfile(repoRoot, 'fun'));
settings = struct;

%% Experimental inputs
settings.inputFile = fullfile(repoRoot, 'input', 'neg', ...
    'checkUp_LGMH1_PAT20261302_G1_Li_T25_Channel_2_Wb_1.xlsx');
settings.capacityFile = fullfile(repoRoot, 'output', 'theoPracCapNeg.mat'); % Q_r [mAh].
settings.cellName = 'LG MH1 PAT half-cell NE';
settings.electrodeName = 'Graphite';
settings.temperatureDegC = 25;
settings.voltageRange = [0 1.5]; % Potential versus Li [V].
settings.activeMaterialFraction = 0.98; % Assumed graphite active-material fraction.
settings.lithiationAnchor = "delithiated"; % Assume theta_min = 0.

%% MSMR fitting inputs
% Six-gallery literature initialization (Verbrugge 2017).
settings.initialMSMR.U0 = [0.08843;0.12799;0.14331;0.16984;0.21446;0.36325]; % V.
settings.initialMSMR.X = [0.43336;0.23963;0.15018;0.05462;0.06744;0.0548]; % Dimensionless.
settings.initialMSMR.w = [0.08611;0.08009;0.72469;2.53277;0.09470;5.97354];
settings.boundScales = [0.5 0.4 0.2]; % Lower-bound factors for U0, X, w.
settings.differentialWeight = 1e-3;
settings.optimizer.seed = 0; % Reproducible GA initialization.
settings.hiddenGalleryCorrection = false;

%% Optional figure settings
% Each inset row is [left bottom width height xmin xmax].
% Positions are normalized; x limits are dimensionless lithiation.
settings.relativeInsets = [.75 .425 .125 .125 .975 1; ...
                          .2 .75 .125 .125 0 .025; ...
                          .45 .5 .125 .125 .45 .55];
settings.fitInsets = [.45 .5 .125 .125 .45 .55; ...
                      .2 .75 .125 .125 .03 .05];
settings.absoluteInsets = [.5 .5 .125 .125 .5 .575; ...
                           .2 .75 .125 .125 .03 .05];
settings.showRelativeLegend = true;
settings.showComparisonLegend = true;

%% Run the shared workflow
settings.outputFile = ""; % Return results without saving a model MAT file.
% To save the model, set an output path before the function call, for example:
% settings.outputFile = fullfile(repoRoot, 'output', 'Uocp_n_example.mat');
[model, ocvrel, supplementaryResults, ocpFigures] = modelOCPelectrode(settings);

% Returned results:
% model                Final electrode model, including MSMR and LUT fields.
% ocvrel               Processed measurements and fitted curves.
% supplementaryResults Snapshot before AM correction.
% ocpFigures           Seven figure handles, indexed in window-number order.
