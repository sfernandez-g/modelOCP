function tests = testGuessMSMR
%TESTGUESSMSMR Gallery ordering must pair each peak with its own capacity.
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
testCase.TestData.originalPath = path;
addpath(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'fun'));
end

function teardownOnce(testCase)
path(testCase.TestData.originalPath);
end

function testAsymmetricGalleriesKeepTheirPeakFractions(testCase)
assumeTrue(testCase, logical(license('test', 'Signal_Toolbox')));
params = struct('U0', [0.15; 0.45], 'X', [0.25; 0.75], 'w', [0.3; 0.3]);
U = linspace(0.05, 0.55, 4001).';
theta = compositionMSMR(U, params, 25);
branch = struct('name', 'charge', 'Z', 1-theta, 'U', U);
other = branch; other.name = 'discharge';
[U0, X, w] = guessMSMR(struct('estimates', [branch, other]), 2, 25);
verifyEqual(testCase, U0, params.U0, 'AbsTol', 0.01);
verifyEqual(testCase, X, params.X, 'AbsTol', 0.03);
verifyTrue(testCase, all(isfinite(w) & w > 0));
end
