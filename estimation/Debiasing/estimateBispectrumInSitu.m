function bisp_est = estimateBispectrumInSitu(shc_samples, L, CGs, f00_est)
% Estimates the bispectrum of signal SHC from noisy, rotated sampled 
% SHC of the signal.
% 
% Input arguments:
%   shc_samples     -   (L+1)^2 by N array of SHC vectors; representing N
%                       samples of SHC vectors of band-limit L
%   L               -   the band-limit of underlying signal
%   CGs             -   Clebsch-Gordan coefficients; 
%                       in the format build_CGs_vec function gives them.
%   f00_est         -   scalar estimating f_{0,0} based on given samples
% 
% Output arguments:
%   bisp_est        -   the estimator of bisp_est
% 
% NOTE:
%   This assumes the noise is in the SHC space, note the image space.
% 

[~, N] = size(shc_samples);

cor_vec = zeros((L+1)^2, 1);
cor_vec(1) = f00_est;
shc_samples = shc_samples - cor_vec;

bisp_size = length(bisp_vec(shc_samples(:, 1), L, CGs));
bisp_est = zeros(bisp_size, 1);
parfor j=1:N
    bisp_est = bisp_est + bisp_vec(shc_samples(:, j), L, CGs);
end

bisp_est = bisp_est/N;
