function bisp_est = estimateBispectrum(bisp_samples, L, f00_est, sigma, DenoVec)
% Estimates the bispectrum based on bispectertum of the samples.
% 
% NOTE:
%   This assumes the noise is in the SHC space, not the image space.
% 

if nargin<5
    DenoVec = buildDenoVec(L);
end

[~, N] = size(bisp_samples);
bisp_est = ( sum(bisp_samples, 2) - sigma^2*conj(f00_est) * DenoVec)/N;
