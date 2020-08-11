function powSpec = powerSpectrum(shc, bandlimit)
%%
% Call format
%   powSpec = powerSpectrum(shc, bandlimit)
% 
% Calculate the power spectrum of the bandlimited spherical function represented by
% spherical harmonics ceofficients in shc
% 
% Input arguments
%   shc             double      (bandlimit+1)^2 x 1 complex array, 
%                               spherical harmonics coefficients.
%   bandlimit       double      positive integer, the bandlimit of shc.
% 
% Output arguments
%   powSpec         double      (bandlimit+1) x 1 array, the power spectrum
%                               of the spherical function represented by
%                               shc, such that
%                                   powSpec(l+1) = l-order spectra of shc.
% 
% Notes
%   (1) This function depends on powerSpectrum_mex.
%   (2) This function performs no input checks.
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
%% Calculate power spectrum
powSpec = powerSpectrum_mex(shc.*conj(shc), bandlimit);
