function b = rotationInvariantBispectrum(coeffs, truncation, angularLimits)
%% 
% Call format
%   b = rotationInvariantBispectrum(coeffs, truncationLimit)
% 
% Compute rotation-invariant bispectrum based on the truncated 
% Fourier-Bessel expansion coefficients coeffs.
% 
% Conventions
%   TODO describe Fourier-Bessel basis convention (as expressed in coeffs)
% 
% Input arguments
%   coeffs          double      array TODO
%   truncation      double      positive scalar, the truncation limit of
%                               the Fourier-Bessel expansion.
%   angularLimits   double      array TODO
% 
% Output arguments
%   b               double      array TODO
% 
% Notes
%   (1) The code performs no input checks.
%   (2) Computing coeffs for this function is best done using the ASPIRE 
%       package implementation of the Fourier-Bessel transform. For
%       information on ASPIRE, see [2].
%   (3) This function depends on rotation_invariant_bispectrum_mex.
%   (4) The general approach for this rotation-invariant bispectrum is
%       described in [1].
% 
% Reference
%   [1] Zhao, Z., & Singer, A. (2014). Rotationally invariant image 
%       representation for viewing direction classification in cryo-EM. 
%       Journal of Structural Biology, 186(1), 153–166. 
%       https://doi.org/https://doi.org/10.1016/j.jsb.2014.03.003
%   [2] http://spr.math.princeton.edu/
% ***********************************************************
% Author    Ido Hadi
% Email     idohadi@mail.tau.ac.il
% Year      2020
% ***********************************************************

%% Compute the rotation-invariant bispectrum
% Complexify the coefficients
% TODO

% Normalize the coefficients
% Explanation: TODO better docs
% See eq. (20) in 
% Rotationally invariant image representation for viewing direction
% classification in cryo-EM
% TODO

% Preprocess the truncation limits
% TODO

% Compute the bispectrum
% TOOD: call a mex function
