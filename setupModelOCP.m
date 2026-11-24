function projectRoot = setupModelOCP()
%SETUPMODELOCP Add this repository's entry points and numerical helpers.
%   Returns the absolute repository path. Call again when switching from a
%   companion project so that this repository's helper versions take priority.
projectRoot = fileparts(mfilename('fullpath'));
addpath(projectRoot);
addpath(fullfile(projectRoot, 'fun'));
end
