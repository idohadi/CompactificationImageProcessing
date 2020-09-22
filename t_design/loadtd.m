function td = loadtd(t)
%%
% Call format
%   td = loadtd(t)
% 
% Load a t-design.
% 
% Input arguments
%   t       double      the t-design.
% 
% Output arguments
%   td      double      n x 3 array, the loaded t-design.
% 
% Notes
%   (1) The t-designs were downloaded from 
%           https://web.maths.unsw.edu.au/~rsw/Sphere/EffSphDes/sf.html
%       Their properties and numerical computation is exlained in [1].
% 
% Reference
%   [1] Womersley, R. S. (2017). Efficient Spherical Designs with Good 
%       Geometric Properties. 
%       https://arxiv.org/abs/1709.01624
% ***********************************************************
% Author    Ido Hadi
% Email     idohadi@mail.tau.ac.il
% Year      2020
% ***********************************************************

dirpath = fullfile(fileparts(which('loadtd.m')), 'tDesigns');
filepath = dir(fullfile(dirpath, ['sf', num2str(t, '%03u'), '.*']));
if isempty(filepath)
    error('No t-design for the given t.');
elseif length(filepath)>1
    error('Multiple t-designs for the given t.');
else
    td = importdata(fullfile(filepath.folder, filepath.name));
end
