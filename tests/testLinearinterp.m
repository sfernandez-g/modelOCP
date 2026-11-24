function tests = testLinearinterp
% Unit tests for linearinterp.
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
testCase.TestData.originalPath = path;
funDir = fullfile(fileparts(fileparts(mfilename('fullpath'))), 'fun');
addpath(funDir);
end

function testinterpolatesAndExtrapolatesUnsortedSamples(testCase)
actual = linearinterp([2; 0; 1], [4; 0; 2], [-1; 0.5; 3]);
verifyEqual(testCase, actual, [-2; 1; 6], 'AbsTol', 1e-12);
end

function testaveragesDuplicateAbscissas(testCase)
actual = linearinterp([0; 1; 1; 2], [0; 1; 3; 4], [1; 1.5]);
verifyEqual(testCase, actual, [2; 3], 'AbsTol', 1e-12);
end

function testpreservesQueryShape(testCase)
actual = linearinterp([0; 1], [0; 2], [0 0.5; 1 1.5]);
verifySize(testCase, actual, [2 2]);
verifyEqual(testCase, actual, [0 1; 2 3], 'AbsTol', 1e-12);
end

function testrejectsInsufficientDistinctSamples(testCase)
verifyError(testCase, @() linearinterp([1; 1], [2; 3], 1), ...
    'linearinterp:InsufficientUniqueData');
end

function testrejectsMismatchedVectors(testCase)
verifyError(testCase, @() linearinterp([0; 1], 1, 0.5), ...
    'linearinterp:InputSizeMismatch');
end

function teardownOnce(testCase)
path(testCase.TestData.originalPath);
end
