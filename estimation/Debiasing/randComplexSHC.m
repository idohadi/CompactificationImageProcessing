function [f, b] = randComplexSHC(L, CGs)
% Generates a random complex SHC vector of band-limit L
% 
% Input arguments
%   L       -   The band-limit of the desired SHC vector
%   CGs     -   Clebsch-Gordan coefficients; 
%               in the format build_CGs_vec function gives them.
% 
% Output arguments:
%   f       -   the SHC vector
%   b       -   its bispectrum
% 

f = rand((L+1)^2, 1);
b = bisp_vec(f, L, CGs);
