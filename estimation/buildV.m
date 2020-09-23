function V = buildV(bandlimit, U)
%%
% Call format
%   V = buildV(L, U)
% 
% Compute the vector V for which
%           1    N  (                               )
%   psEst = - * sum ( powerSpectrum(n-th sample SHC) - sigma^2 * V )
%           N   n=1 (                               )
% is an unbiased estimator.
% 
% Input arguments
%   bandlimit   double      positive integer, the bandlimit of the sampled
%                           SHCs.
%   U           double      matrix, the output of buildU. See description
%                           there.
% 
% Output arguments
%   V               double      (bandlimit+1) x 1 array, the vector
%                               described above.
% 
% Notes
%   This function performs no input checks.
% 
% Reference
%   None
% ***********************************************************
% Author    Ido Hadi
% Email     idohadi@mail.tau.ac.il
% Year      2020
% ***********************************************************

%% Compute V
V = buildV_mex(U*U', bandlimit);
V = real(V);
