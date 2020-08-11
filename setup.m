%% Setup file
[p1, ~, ~] = fileparts(mfilename('fullpath'));
cd(p1);

%% Run external setup files
run(fullfile(p1, 'extern', 'FastSphericalHarmonicsTransform', 'setup.m'));
run(fullfile(p1, 'extern', 'SmallRotationToolbox', 'setup.m'));

%% Adds the necessary folders to path
addpath(p1);
addpath(genpath(fullfile(p1, 'alignment')));
addpath(genpath(fullfile(p1, 'clebsch_gordan')));
addpath(genpath(fullfile(p1, 'spherical_harmonics')));
addpath(genpath(fullfile(p1, 'spherical_spectra')));