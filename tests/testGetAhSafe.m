function tests = testGetAhSafe
tests = functiontests(localfunctions);
end

function setupOnce(~)
funDir = fullfile(fileparts(fileparts(mfilename('fullpath'))), 'fun');
addpath(funDir);
end

function testshiftsSelectedCapacityToZero(testCase)
actual = getAhSafe([0; 0.1; 0.3; 0.8], [2; 3; 4]);
verifyEqual(testCase, actual, [0; 0.2; 0.7], 'AbsTol', 1e-12);
end

function testemptySelectionReturnsZero(testCase)
verifyEqual(testCase, getAhSafe([0; 1], []), 0);
end
