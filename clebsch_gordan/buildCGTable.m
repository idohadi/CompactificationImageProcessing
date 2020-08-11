function CGTable = buildCGTable(bandlimit, filename)
%%
% Call format
%   CGTable = buildCGTable(bandlimit)
%   CGTable = buildCGTable(bandlimit, filename)
% 
% Calculate all Clebsch-Gordan coefficients for a given bandlimit and save
% them.
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
% 
% Optional arguments
%   filename        char        path to file relative to currently open
%                               folder.
% 
% Notes
%   (1) This function depends on buildCGTable_mex.
%   (2) The coefficients are calcualted using  the method of [1]. In [1]
%       there is also a concise description of theoretical background.
%   (3) It is possible to be more memory efficient, by using symmetries
%       inside the Clebsch-Gordan coefficients.
%   (4) If filename is specified, it is best to use an absolute path.
% 
% Reference
%   [1] Straub, W. O. (n.d.). Efficient Computation of Clebsch-Gordan 
%       Coefficients. Retrieved October 28, 2019, 
%       from http://vixra.org/abs/1403.0263
% ***********************************************************
% Author    Ido Hadi
% Email     idohadi@mail.tau.ac.il
% Year      2020
% ***********************************************************
%% Input validation
assert(round(bandlimit)==bandlimit && bandlimit>=1, ...
    'The bandlimit must be a positive integer.');

if nargin==1
    filename = ['ClebschGordanCoeffs/CGT', num2str(bandlimit), '.mat'];
end

%% Build the table
CGTable = buildCGTable_mex(bandlimit);

%% Save it
saveCGTable(CGTable, filename);
