function shc = randSHC(bandlimit, N)
%%
% Call format
%   shc = randSHC(bandlimit)
%   shc = randSHC(bandlimit, N)
% 
% Uniformly sample the spherical harmonics coefficients having power
% spectrum of all ones and representing a function with bandlimit
% bandlimit.
% 
% Input arguments
%   bandlimit       double      positive integer, the bandlimit of shc.
% 
% Output arguments
%   nshc            double      (2*N)^2 x 1 complex array, spherical 
%                               harmonics coefficients, normalized such
%                               that it has a power spectrm of all ones.
% 
% Optional arguments
%   Name            Default         Definition
%   N               1               Number of spherical harmonics
%                                   coefficients to sample.
% 
% Notes
%   (1) This function performs no input checks.
%   (2) This function depends on normSHC_mex.
% 
% Reference
%   None
% ***********************************************************
% Author    Ido Hadi
% Email     idohadi@mail.tau.ac.il
% Year      2020
% ***********************************************************
%% Input validation
if nargin<2
    N = 1;
end

%% Sampling SHCs
shc = randn((bandlimit+1)^2, N);
shc = normSHC(shc, bandlimit);
