function loadCGTable(bandlimit)
%%
% Call format
%   loadCGTable(bandlimit)
% 
% Load a saved Clebsch-Gordan table from a .mat file. The loaded table is
% saved in a global variable called CGs.
% 
% CGTable format specification
%   For 
%       0<=l1<=bandlimit, 
%       0<=l2<=l2, 
%       abs(l1-l2)<=l<=min(l1+l2,bandlimit),  and
%       abs(m)<=l.
%   
%   The table has the following format:
%   Name                                                                    Type        Size        
%   CGTable                                                                 cell        (bandlimit+1) x 1
%   CGTable{l1+1}                                                           cell        (l1+1) x 1
%   CGTable{l1+1}{l2+1}                                                     cell        (min(l1+l2, bandlimit)-abs(l1-l2)+1) x 1
%   CGTable{l1+1}{l2+1}{l1+l2-abs(l1-l2)+1}                                 cell        (2*l+1) x 1
%   CGTable{l1+1}{l2+1}{l1+l2-abs(l1-l2)+1}{m+l+1}                          double      (max{-l1, m-l2} - min{l1, m+l2} + 1) x 1
%   CGTable{l1+1}{l2+1}{l1+l2-abs(l1-l2)+1}{m+l+1}(m1+max{-l1, m-l2}+1)     double      1 x 1
%   
%   Overall,
%       CGTable{l1+1}{l2+1}{l1+l2- abs(l1-l2)+1}{m+l+1}(m1+max{-l1, m-l2}+1)
%           = < l1 m1 l2 (m-m1) | l m >
% 
% Input arguments
%   bandlimit       double      self-explanatory.
% 
% Output arguments
%   None
% 
% Notes
%   None
% 
% Reference
%   None
% ***********************************************************
% Author    Ido Hadi
% Email     idohadi@mail.tau.ac.il
% Year      2020
% ***********************************************************

%% Input validation
assert(round(bandlimit)==bandlimit && bandlimit>=1, ...
    'The bandlimit must be a positive integer.');

%% Load CG file
[rfn, ~, ~] = fileparts(mfilename('fullpath'));
filename = fullfile(rfn, 'ClebschGordanCoeffs', ['CGT', num2str(bandlimit), '.mat']);

global CGs
CGs = load(filename, '-mat');
CGs = CGs.CGTable;
