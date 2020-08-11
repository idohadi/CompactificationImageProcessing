function nshc = normSHC(shc, bandlimit)
%%
% Call format
%   nshc = normSHC(shc, bandlimit)
% 
% Normalize spherical harmonics coefficients array to have a power spectrum
% of all ones.
% 
% Input arguments
%   shc             double      (bandlimit+1)^2 x N complex array, 
%                               spherical harmonics coefficients.
%   bandlimit       double      positive integer, the bandlimit of shc.
% 
% Output arguments
%   nshc            double      (2*N)^2 x 1 complex array, spherical 
%                               harmonics coefficients, normalized such
%                               that it has a power spectrm of all ones.
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
%% Normalizing spherical harmonics coefficients
nshc = normSHC_mex(shc, real(shc.*conj(shc)), bandlimit);
