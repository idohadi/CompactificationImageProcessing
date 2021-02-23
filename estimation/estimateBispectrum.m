function bispEst = esimateBispectrum(batched, images, bandlimit, tDesign, interval, scalingParam, debiasingMatrix, sigma2)
%%
% Call format
%   bispEst = esimateBispectrum(images, bandlimit, tDesign, interval, scalingParam, debiasingMatrix, sigma2)
% 
% Estimate the bispectrum from images by
%           1    N  (                                                                         )
% bispEst = - * sum ( bispectrum(n-th sample SHC) - sigma^2 * K * (realified n-th sample SHC) )
%           N   n=1 (                                                                         )
% where K = debiasingMatrix.
% 
% If batched is true, parfor will be used.Otherwise, a regular for loop
% will be used.
% 
% Input arguments
%   batched             logical     true if parfor will be used to shared
%                                   the workload among serveral workers.
%                                   false if a single worker will be
%                                   utilized.
%   images              double      imageSize x imageSize x sampleSize
%                                   array, the images dataset.
%   bandlimit           double      positive integer, the bandlimit of the 
%                                   sampled SHCs.
%   tDesign             double      N x 3 array, a spherical design (a 
%                                   t-design) in Cartesian coordiantes.
%   interval            double      1 x 2 array, interval(1) is the lower 
%                                   bound of the interval and interval(2).
%   scalingParam        double      positive number, scaling parameter for 
%                                   the projection.
%   debiasingMatrix     double      imageSize^2 x 2*(bandlimit+1)^2 array, 
%                                   the debiasing matrix.
%   sigma2              double      positive scalar, the noise variance, 
%                                   sigma^2. 
% 
% Output arguments
%   bispEst             double      k x 1 array, the estimator of the 
%                                   bispectrum.
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

%% Estimate the bispectrum
% TODO
