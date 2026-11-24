function tests = testFitMSMR
%TESTFITMSMR Small synthetic optimizer check, independent of private XLSX data.
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
testCase.TestData.originalPath = path;
addpath(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'fun'));
end

function teardownOnce(testCase)
path(testCase.TestData.originalPath);
end

function testSyntheticFitAndRandomStream(testCase)
assumeTrue(testCase, license('test', 'Optimization_Toolbox') ...
    && license('test', 'GADS_Toolbox'));
params = struct('U0', 0.25, 'X', 1, 'w', 0.5, 'thetamin', 0.1, 'thetamax', 0.9);
U = linspace(0.20, 0.30, 41).';
[theta, derivative] = compositionMSMR(U, params, 25);
estimate = struct('name', 'charge', 'U', U, 'Z', 1-(theta-0.1)/0.8, 'dZ', derivative/0.8);
% Use measured endpoints to define the interval exactly as the fitter does.
params.thetamin = min(theta); params.thetamax = max(theta);
config = struct('init', params, 'lb', params, 'ub', params, 'optw', 1e-2);
config.optimizer = struct('seed', 2027, ...
    'gaOptions', {{'Display', 'off', 'PopulationSize', 10, 'MaxGenerations', 2}}, ...
    'fminconOptions', {{'Display', 'off'}});
before = rng;
[fitted, composition, ~, ~, ~, ~, ~, inverse, info] = ...
    fitMSMR(config, struct('estimates', estimate), 25, "charge", 0, struct('A', [], 'b', []));
verifyEqual(testCase, rng, before);
verifyEqual(testCase, fitted, params);
verifyEqual(testCase, composition, theta, 'AbsTol', 1e-12);
verifyEqual(testCase, inverse(theta(21)), U(21), 'AbsTol', 1e-8);
verifyLessThanOrEqual(testCase, info.maxConstraintViolation, 1e-10);
end
