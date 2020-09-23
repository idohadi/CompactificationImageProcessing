function shc = randSHCReal(bandlimit, N)
%%
% Call format
%   shc = randSHCReal(bandlimit)
%   shc = randSHCReal(bandlimit, N)
% 
% Uniformly sample the spherical harmonics coefficients having power
% spectrum of all ones and representing a real-valued function with 
% bandlimit bandlimit.
% 
% Input arguments
%   bandlimit       double      positive integer, the bandlimit of shc.
% 
% Output arguments
%   shc             double      (bandlimit+1)^2 x N complex array, spherical 
%                               harmonics coefficients, normalized such
%                               that it has a power spectrm of all ones, 
%                               representing a real-valued, bandlimited
%                               spherical function.
% 
% Optional arguments
%   Name            Default         Definition
%   N               1               Number of spherical functions to sample.
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

%% Input validation
if nargin<2
    N = 1;
end

%% Sampling SHCs
tmp = randn((bandlimit+1)^2, N);
shc = complex(zeros((bandlimit+1)^2, 1));

% Handle the spherical harmonic coefficient of order 0 and degre 0
shc(1) = tmp(1);

% Handle the spherical harmonic coefficients of order greater than 0
for l=1:bandlimit
    % Handle negative degress
    
    shc(l*(l+1)+(-l:-1)+1) = ((-1).^(-l:-1)').*flip((tmp(l*(l+1)+(-l+1:2:l)+1) ...
                            - 1i*tmp(l*(l+1)+(-l+2:2:l)+1)));
    
    % Handle degree zero
    shc(l*(l+1)+1) = tmp(l^2+1);
    
    % Handle positive degress
    shc(l*(l+1)+(1:l)+1) = tmp(l*(l+1)+(-l+1:2:l)+1) ...
                            + 1i*tmp(l*(l+1)+(-l+2:2:l)+1);
end

shc = normSHC(shc, bandlimit);
