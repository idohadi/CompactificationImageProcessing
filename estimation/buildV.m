function V = buildV(bandlimit, U)
%%
% Call format
%   V = buildV(L, U)
% 
% Compute the vector V for which
%           1    N  (                               )
%   psEst = - * sum ( powerSpectrum(n-th sample SHC) - V )
%           N   n=1 (                               )
% is an unbiased estimator.
% 
% Input arguments
%   U           double      
%   bandlimit   double      positive integer, the bandlimit of the sampled
%                           SHCs.
% 
% Output arguments
%   P               double      a matrix represneting the linear opeartor
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
