function tests = testCompositionMSMR
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
testCase.TestData.originalPath = path;
funDir = fullfile(fileparts(fileparts(mfilename('fullpath'))), 'fun');
addpath(funDir);
end

function teardownOnce(testCase)
path(testCase.TestData.originalPath);
end

function testcompositionIsBoundedAndDecreasesWithPotential(testCase)
params.U0 = [0.1; 0.3; 0.5];
params.X = [0.2; 0.5; 0.3];
params.w = [0.08; 0.12; 0.1];
U = linspace(-0.2, 0.8, 201)';
[Z, dZ] = compositionMSMR(U, params, 25);
verifyTrue(testCase, all(Z >= 0 & Z <= sum(params.X)));
deltaZ = diff(Z);
verifyTrue(testCase, all(deltaZ <= 1e-14) && any(deltaZ < -1e-14));
verifyTrue(testCase, all(isfinite(dZ)));
end

function testdifferentialMatchesFiniteDifferenceMagnitude(testCase)
params.U0 = [0.15; 0.4];
params.X = [0.35; 0.65];
params.w = [0.1; 0.16];
U = linspace(-0.1, 0.7, 1001)';
[Z, dZ] = compositionMSMR(U, params, 25);
fd = abs(diff(Z) ./ diff(U));
verifyEqual(testCase, dZ(2:end-1), (fd(1:end-1) + fd(2:end))/2, ...
    'RelTol', 0.025, 'AbsTol', 1e-5);
end

function testSaturatedGalleriesRemainFinite(testCase)
params = struct('U0', [0.15; 0.4], 'X', [0.35; 0.65], 'w', [0.1; 0.16]);
U = [-1e4, 1e4];
[Z, dZ, xj, Uj, dZj] = compositionMSMR(U, params, 25);
verifyEqual(testCase, Z, [1; 0], 'AbsTol', 1e-14);
verifyEqual(testCase, dZ, [0; 0], 'AbsTol', 1e-14);
verifyEqual(testCase, Uj, repmat(U, 2, 1));
verifyTrue(testCase, all(isfinite([xj(:); dZj(:); Uj(:)])));
end

function testGalleryDerivativesSumToTotal(testCase)
params = struct('U0', [0.15; 0.4], 'X', [0.35; 0.65], 'w', [0.1; 0.16]);
[~, dZ, ~, ~, dZj] = compositionMSMR(linspace(0, 0.6, 501), params, 25);
f = 96485.3329/(8.3144598*(25+273.15));
verifyEqual(testCase, dZ, f*sum(dZj, 1).', 'AbsTol', 1e-12);
end
