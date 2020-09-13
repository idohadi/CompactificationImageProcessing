function f00_est = estimatef00(shc_samples)
% Estimate the f_{0,0} given samples of SHC.
% 
% Input parameters:
%   shc_samples     -   (L+1)^2 by N array of SHC vectors; representing N
%                       samples of SHC vectors of band-limit L
% 
% Output parameters:
%   f00_est         -   scalar estimating f_{0,0} based on given samples
% 
% NOTE:
%   This assumes the noise is in the SHC space, not the image space.
% 

[~, N] = size(shc_samples);
f00_est = sum(shc_samples(1, :))/N;
