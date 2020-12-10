function b = rotationInvariantBispectrum(coeff, truncation, angularLimits, bLen)
%% 
% Call format
%   b = rotationInvariantBispectrum(coeff, truncationLimit)
% 
% Compute rotation-invariant bispectrum based on the truncated 
% Fourier-Bessel expansion coefficients coeff.
% 
% Conventions
%   TODO describe Fourier-Bessel basis convention (as expressed in coeff)
% 
% Input arguments
%   coeff           double      array TODO
%   truncation      double      positive scalar, the truncation limit of
%                               the Fourier-Bessel expansion.
%   angularLimits   double      array TODO
% 
% Output arguments
%   b               double      array TODO
% 
% Notes
%   (1) The code performs no input checks.
%   (2) Computing coeff for this function is best done using the ASPIRE 
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
angularLimits = angularLimits(:);

tmp = angularLimits(1) + [0; cumsum(2*angularLimits(2:end))];

tmp2 = cumsum(angularLimits(2:end));
tmp2 = [0; tmp2];
coeffComp = zeros(sum(angularLimits(2:end)), 1);

for J=1:truncation
    c = coeff(tmp(J)+1:tmp(J+1));
    c = reshape(c, [angularLimits(J+1), 2]);
    coeffComp(tmp2(J)+1:tmp2(J+1)) = 0.5*(c(:, 1) + 1i*c(:, 2));
end

% Normalize the coefficients
% Explanation: TODO better docs
% See eq. (20) in 
% Rotationally invariant image representation for viewing direction
% classification in cryo-EM
coeffComp = coeffComp./(abs(coeffComp).^(2/3));

% Compute the bispectrum
b = zeros(bLen + angularLimits(1), 1);
b(1:angularLimits(1)) = coeff(1:angularLimits(1));

m = angularLimits(1)+1;
for k1=1:truncation
    for k2=1:min(k1, truncation-k1)
        M = coeffComp(tmp2(k1)+1:tmp2(k1+1)) * coeffComp(tmp2(k2)+1:tmp2(k2+1)).';
        M = M(:);
        
        c = conj(coeffComp(tmp2(k1+k2)+1:tmp2(k1+k2+1)));        
        c = c.';
        
        M = M .* c;
        M = M(:);
        
        b(m:m+length(M)-1) = M;
        m = m + length(M);
    end
end
