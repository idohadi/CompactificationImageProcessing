function K = buildK(bandlimit, U)
%%
% Call format
%   K = buildK(bandlimit, U)
% 
% Compute the matrix V for which
%           1    N  (                                                                         )
%   psEst = - * sum ( bispectrum(n-th sample SHC) - sigma^2 * K * (realified n-th sample SHC) )
%           N   n=1 (                                                                         )
% is an unbiased estimator.
% 
% Input arguments
%   bandlimit   double      positive integer, the bandlimit of the sampled
%                           SHCs.
%   U           double      matrix, the output of buildU. See description
%                           there.
% 
% Output arguments
%   K               double      bispectrum legnth x (2*(bandlimit+1)^2) 
%                               array, the matrix described above.
% 
% Notes
%   (1) This function performs no input checks.
%   (2) This function depends on buildK_mex.
% 
% Reference
%   None
% ***********************************************************
% Author    Ido Hadi
% Email     idohadi@mail.tau.ac.il
% Year      2020
% ***********************************************************

%% Compute K
UUH = full(U*U');
UUT = full(U*U.');

K = buildK_mex(UUH, UUT, bandlimit);
K = sparse(K.');
