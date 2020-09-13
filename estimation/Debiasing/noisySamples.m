function shc_samples = noisySamples(f, L, N, sigma, rotate)
% Creates N samples of noisy SHC f of band-limit L. 
% The noise is Gaussian, with zero mean and variance sigma.
% 
% Input arguments:
%   f               -   the original SHC; the ground truth.
%   L               -   the band-limit of f
%   N               -   The numbero of samples to generate
%   sigma           -   The variance of the Gaussian noise
%   rotate          -   boolian; 
%                       if true, rotate the signal using uniformly chosen
%                       Euler angles. 
%                       If false, do not rotate.
% 
% Output arguments:
%   shc_samples     -   (L+1)^2 by N array of samples
% 

if nargin<=4
    rotate = false;
end

if ~rotate
    shc_samples = f + sigma*(randn((L+1)^2, N) + 1i*randn((L+1)^2, N))/sqrt(2);
else
%     Create random rotation angles
    alpha = 2*pi*rand(N, 1);
    beta = 2*pi*rand(N, 1);
    gamma = 2*pi*rand(N, 1);
    
%     Apply rotations to original SHC f
    shc_samples = zeros((L+1)^2, N);
    for j=1:N
        shc_samples(:, j) = rotshc(f, alpha(j), beta(j), gamma(j));
    end
    
%     Add noise
    shc_samples = shc_samples ...
        + sigma*(randn((L+1)^2, N) + 1i*randn((L+1)^2, N))/sqrt(2);
end
