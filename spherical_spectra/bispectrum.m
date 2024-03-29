function [b, grad] = bispectrum(shc, bandlimit, CGs)
%%
% Call format
%   b = bispectrum(shc, bandlimit)
%   [b, grad] = bispectrum(shc, bandlimit)
% 
% Calculate the bispectrum and the gradient of the spherical harmonice
% coefficients.
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
%   shc             double      (bandlimit+1)^2 x 1 complex array, 
%                               spherical harmonics coefficients.
%   bandlimit       double      positive integer, the bandlimit of shc.
% 
% Output arguments
%   b               double      bispectrum f spherical function
%                               represented by shc.
%   grad            double      gradient of b as a function of 
%                                   [real(shc(1)); imag(shc(1)); ...
%                                   real(shc(end)); imag(shc(end))]
% 
% Notes
%   (1) This function depends on bispectrum_mex.
%   (2) This function performs no input checks.
% 
% Reference
%   None
% ***********************************************************
% Author    Ido Hadi
% Email     idohadi@mail.tau.ac.il
% Year      2020
% ***********************************************************

%% Calculate the bispectrum
% global CGs;
if nargout==1
    b = bispectrum_mex(shc, bandlimit, CGs);
elseif nargout==2
    [b, gradI, gradJ, gradVals, gradNNZ, gradRowsNo, gradColsNo] ...
        = bispectrum_mex(shc, bandlimit, CGs);
    gradVals(isnan(gradVals)) = 0;
    grad = sparse(gradI, gradJ, gradVals, gradRowsNo, gradColsNo, gradNNZ);
end
b(isnan(b)) = 0;
