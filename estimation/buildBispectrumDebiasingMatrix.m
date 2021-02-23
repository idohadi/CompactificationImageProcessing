function debiasingMatrix = buildBispectrumDebiasingMatrix(imageSize, bandlimit, tDesign, interval, scalingParam)
%%
% Call format
%   debiasingMatrix = buildBispectrumDebiasingMatrix(imageSize, tDesign, interval, scalingParam)
% 
% Build the debiasing matrix for bispectrum estimation from images of size
% imageSize x imageSize. This is a matrix K for which
%           1    N  (                                                                         )
% bispEst = - * sum ( bispectrum(n-th sample SHC) - sigma^2 * K * (realified n-th sample SHC) )
%           N   n=1 (                                                                         )
% is an unbiased estimator.
% 
% Input arguments
%   imageSize       double      positive integer, the bispectrum is
%                               calcualted for imageSize x imageSize 
%                               images.
%   bandlimit       double      positive integer, the bandlimit of the sampled
%                               SHCs.
%   tDesign         double      N x 3 array, a spherical design (a 
%                               t-design) in Cartesian coordiantes.
%   interval        double      1 x 2 array, interval(1) is the lower bound
%                               of the interval and interval(2).
%   scalingParam    double      positive number, scaling parameter for 
%                               the projection.
% 
% Output arguments
%   debiasingMatrix double      imageSize^2 x 2*(bandlimit+1)^2 array, 
%                               the debiasing matrix.
% 
% Notes
%   None
% 
% Reference
%   None
% ***********************************************************
% Author    Ido Hadi
% Email     idohadi@mail.tau.ac.il
% Year      2021
% ***********************************************************

%% Build the debiasing matrix
U = buildU(bandlimit, imageSize, tDesign, interval, scalingParam);
debiasingMatrix = buildK(bandlimit, U);
