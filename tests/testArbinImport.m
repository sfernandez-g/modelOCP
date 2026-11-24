function tests = testArbinImport
%TESTARBINIMPORT Check integration/cache behavior using tiny synthetic workbooks.
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
testCase.TestData.originalPath = path;
root = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root, 'fun'));
buildDirectory = fullfile(root, '.build');
if ~isfolder(buildDirectory), mkdir(buildDirectory); end
testCase.TestData.directory = tempname(buildDirectory);
mkdir(testCase.TestData.directory);
end

function teardownOnce(testCase)
path(testCase.TestData.originalPath);
% Retain tiny synthetic fixtures under the ignored .build directory for audit.
% Cloud-synced Windows folders can keep workbook directories temporarily locked.
end

function testCapacityAndSeparateCache(testCase)
file = fullfile(testCase.TestData.directory, 'synthetic_Wb_1.xlsx');
writeFixture(file, [0; 1; 2; 3]);
data = getDataArbin_XLSX(file, 7);
verifyEqual(testCase, data.ch_capacity(end), 3/3600, 'AbsTol', 1e-14);
verifyEqual(testCase, data.dch_capacity, zeros(4, 1));
verifyTrue(testCase, isfile(fullfile(testCase.TestData.directory, '.cache', 'synthetic_Wb_1.mat')));
verifyFalse(testCase, isfile(fullfile(testCase.TestData.directory, 'synthetic_Wb_1.mat')));
cached = getDataArbin_XLSX(file, 7);
verifyEqual(testCase, cached.capacity, data.capacity);
end

function testTimeResetIsRejected(testCase)
file = fullfile(testCase.TestData.directory, 'reset_Wb_1.xlsx');
writeFixture(file, [0; 2; 1]);
verifyError(testCase, @() getDataArbin_XLSX(file, 7), 'getDataArbin_XLSX:TimeReset');
end

function testMissingSourceIsReported(testCase)
verifyError(testCase, @() getDataArbin_XLSX(fullfile(testCase.TestData.directory, 'absent.xlsx')), ...
    'getDataArbin_XLSX:MissingSource');
end

function writeFixture(file, time)
metadata = array2table(ones(1, 10), 'VariableNames', ...
    {'Channel', 'Start', 'A', 'Fixture', 'Test Name', 'B', 'C', 'Var8', 'D', 'Var10'});
writetable(metadata, file, 'Sheet', 'Info');
measurements = table(time, ones(size(time)), 3.5*ones(size(time)), 7*ones(size(time)), ...
    'VariableNames', {'Test Time (s)', 'Current (A)', 'Voltage (V)', 'Step Index'});
writetable(measurements, file, 'Sheet', 'Measurements');
end
