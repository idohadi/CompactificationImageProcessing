function saveCGTable(CGTable, filename)
%%
% save the CG table to file
% Default is to save a .mat file
% Alternatively, save a csv file with a table with headers
%       l1 m1 l2 m2 | l m   value
% 
%%
% Call format
%   saveCGTable(CGTable, filename)
% 
% Save Clebsch-Gordan table to file as a .mat file.
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
%   CGTable         double      A cell array. Its specification is above.
%   filename        char        path to file relative to currently open
%                               folder.
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

%% Save table
save(filename, 'CGTable', '-mat');
