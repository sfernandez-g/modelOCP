function tests = testProcessCU
%TESTPROCESSCU Synthetic phases check orientation, units and branch labels.
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
testCase.TestData.originalPath = path;
addpath(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'fun'));
end

function teardownOnce(testCase)
path(testCase.TestData.originalPath);
end

function testRowsAndColumnsGiveSameCapacity(testCase)
[data, cfg] = fixture();
expected = processCU(data, cfg, {'MainVolt'}, 1);
data.step = data.step.'; % A mixed row/column mask must not expand to N-by-N.
actual = processCU(data, cfg, {'MainVolt'}, 1);
verifyEqual(testCase, actual.QAh, 1, 'AbsTol', 1e-12);
verifyEqual(testCase, actual, expected);
end

function testTemperatureAndEstimatorNames(testCase)
assumeTrue(testCase, logical(license('test', 'Signal_Toolbox')));
[data, cfg] = fixture();
[cal, ocv, names] = processOCPCurves(data, cfg, {'MainVolt'}, 1, 'synthetic', ...
    [2.5, 4.2], @(z) z > 0.1 & z < 0.9, 501, 30);
verifyEqual(testCase, cal.TdegC, 30);
verifyEqual(testCase, ocv.TdegC, 30);
verifyEqual(testCase, string(names), ["discharge", "charge", "diagonalAverage"]);
verifyEqual(testCase, string({ocv.estimates.name}), string(names));
end

function [data, cfg] = fixture()
n = 101;
data = struct('step', repelem([5; 6; 3; 4], n), 'TC_Counter1', ones(4*n, 1), ...
    't', (0:4*n-1).', 'i', repelem([-1; -0.1; 1; 0.1], n), ...
    'dch_capacity', [linspace(0, 0.8, n).'; linspace(0.8, 1, n).'; ones(2*n, 1)], ...
    'ch_capacity', [zeros(2*n, 1); linspace(0, 0.8, n).'; linspace(0.8, 1, n).'], ...
    'v', [linspace(4.2, 2.5, n).'; 2.5*ones(n, 1); linspace(2.5, 4.2, n).'; 4.2*ones(n, 1)]);
cfg = struct('name', {'CC_dis', 'CV_dis', 'CC_chg', 'CV_chg'}, ...
    'step', {5, 6, 3, 4}, 'counter', {1, 1, 1, 1});
end
