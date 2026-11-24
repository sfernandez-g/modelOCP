% PROCESSOCPCURVES Run the shared capacity, calibration, and OCV estimation steps.
%
% [cal,ocvrel,estimators] = processOCPCurves(testData,cfg,channels,channelIndex,...
%     cellName,voltageFilter,socFilter,resolution)
%
% Plotting, absolute-lithiation assumptions, and model fitting remain the
% responsibility of the calling workflow.
function [cal, ocvrel, estimators] = processOCPCurves(testData, cfg, channels, channelIndex, cellName, voltageFilter, socFilter, resolution, temperatureDegC)

% The optional ninth argument is temperature [degC], default 25.
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
