function results = runModelOCPTests()
%RUNMODELOCPTESTS Run unit and small synthetic-fit tests without private data.
%   Solver tests are skipped when Optimization/Global Optimization are absent.
originalPath = path;
cleanup = onCleanup(@() path(originalPath));
root = setupModelOCP();
results = runtests(fullfile(root, 'tests'));
assertSuccess(results);
end
