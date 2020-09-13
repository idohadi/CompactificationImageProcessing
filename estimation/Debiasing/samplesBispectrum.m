function bisp_samples = samplesBispectrum(shc_samples, L, CGs, f00_est)
% Calculate the bispectrum of f - the estimate of f_{0,0}.
% 
% Input arguments:
%   shc_samples     -   (L+1)^2 by N array of SHC vectors; representing N
%                       samples of SHC vectors of band-limit L
%   L               -   The band-limit of the samples; integer
%   CGs             -   Clebsch-Gordan coefficients; 
%                       in the format build_CGs_vec function gives them.
%   f00_est         -   scalar estimating f_{0,0} based on given samples
% 
% Output argument:
%   bisp_samples    -   N columns array; the bispectrum of 
%                       each of the samples
% 

[~, N] = size(shc_samples);

if nargin>3
    cor_vec = zeros((L+1)^2, 1);
    cor_vec(1) = f00_est;
    shc_samples = shc_samples - cor_vec;
end

% Bispectrum of sample SHCs
bisp_samples = cell(N, 1); 
parfor j=1:N
    bisp_samples{j} = bisp_vec(shc_samples(:, j), L, CGs);
end
bisp_samples = reshape(cell2mat(bisp_samples), length(bisp_samples{1}), N);
