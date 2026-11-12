function tests = testCompositionMSMR
tests = functiontests(localfunctions);
end

function setupOnce(~)
funDir = fullfile(fileparts(fileparts(mfilename('fullpath'))), 'fun');
addpath(funDir);
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
