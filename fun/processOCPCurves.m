% PROCESSOCPCURVES Run the shared capacity, calibration, and OCV estimation steps.
%
% [cal,ocvrel,estimators] = processOCPCurves(testData,cfg,channels,channelIndex,...
%     cellName,voltageFilter,socFilter,resolution)
% [cal,ocvrel,estimators] = processOCPCurves(testData,cfg,channels,channelIndex,...
%     cellName,voltageFilter,socFilter,resolution,temperatureDegC)
%
% testData is the Arbin import structure. cfg, channels, and channelIndex
% select the test steps and voltage channel as documented in processCU.
% cellName identifies the electrode. voltageFilter is [min max] [V];
% socFilter maps a dimensionless SOC vector to a logical selection mask.
% Both filters select the intervals used for diagonal-average correlation.
% resolution is the number of interpolation points. temperatureDegC is
% optional and defaults to 25 [degC].
%
% cal contains the processed branches and reference capacity QAh [Ah].
% ocvrel contains the discharge, charge, and diagonal-average estimates.
% estimators is a cell array of their names in the returned order.
%
% Plotting, absolute-lithiation assumptions, and model fitting remain the
% responsibility of the calling workflow.
function [cal, ocvrel, estimators] = processOCPCurves(testData, cfg, channels, channelIndex, cellName, voltageFilter, socFilter, resolution, temperatureDegC)

if nargin < 9, temperatureDegC = 25; end
validateattributes(temperatureDegC, {'numeric'}, {'scalar', 'real', 'finite', '>', -273.15});
cal = processCU(testData, cfg, channels, channelIndex);
cal.TdegC = temperatureDegC;
cellData.raw = cal;
cellData.cal = cal;
cellData.name = cellName;

config.vmin = cal.minV;
config.vmax = cal.maxV;
config.datype = 'useAvg';
config.daxcorrfiltz = @(Z) {socFilter(Z)};
config.daxcorrfiltv = @(V) {voltageFilter(1) < V & V < voltageFilter(2)};

estimators = {'diagonalAverage'};
ocvrel = estimateOCV(cellData, config, resolution, estimators);
ocvrel = addDchChOCV(ocvrel, cellData, "charge"); estimators{end+1} = 'charge';
ocvrel = addDchChOCV(ocvrel, cellData, "discharge"); estimators{end+1} = 'discharge';
ocvrel.Q = cal.QAh;
ocvrel.estimates = ocvrel.estimates([3, 2, 1]);
% Keep names aligned with the returned data; callers use these names to index LUTs.
estimators = {ocvrel.estimates.name};
end
