function bispEst = estimateBispectrum(batched, images, bandlimit, tDesign, interval, scalingParam, debiasingMatrix, sigma2, sh)
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
%   sh                  double      spherical harmonics matrix.
% 
% Output arguments
%   bispEst             double      k x 1 array, the estimator of the 
%                                   bispectrum.
% 
% Notes
%   (1) Prior to using this function, the Clebsch-Gordan coefficient table
%       needs to be loaded as a global variable. Use the loadCGTable
%       function for that.
% 
% Reference
%   None
% ***********************************************************
% Author    Ido Hadi
% Email     idohadi@mail.tau.ac.il
% Year      2021
% ***********************************************************

%% Input validation
% batched
assert(isscalar(batched) && islogical(batched), ...
    'batched must be a logical scalar.');
assert(size(images, 1)==size(images, 2), ...
    'Images must be rectangular.');
assert(isscalar(bandlimit) && bandlimit>0 && round(bandlimit)==bandlimit, ...
    'Bandlimit must be a positive integer.');
assert(size(tDesign, 2)==3, ...
    't-design must be a 3 x N array.');
assert(size(interval, 1)==1 && size(interval, 2), ...
    'Interval must be a 1 x 2 array.');
assert(isscalar(scalingParam) && scalingParam>=0, ...
    'Scaling parameter must be a non-negative scalar.');
assert(isscalar(sigma2) && sigma2>=0, ...
    'Variance must be a non-negative scalar.');

sampleSize = size(images, 3);
imageSize = size(images, 1);

%% Estimate the bispectrum
% Declare global CGs
global CGs;

% Get bispectrum vector length
shc = image2shc(images(:, :, 1), bandlimit, tDesign, interval, scalingParam);
b = bispectrum(shc, bandlimit, CGs);
bispLen = size(b, 1);
clear shc b;

if batched
    % Turn on parallel pool
    gcp;
    
    % Set up batching
    batchSize = 10^3;
    batchNo = ceil(sampleSize/batchSize);
    
    batchSampleSize = repmat(batchSize, 1, batchNo);
    if batchSize*batchNo > sampleSize
        batchSampleSize(end) = mod(sampleSize, batchNo);
    end
    
    batchUpperBound = cumsum(batchSampleSize);
    batchLowerBounds = [1, batchUpperBound(1:end-1)];
    
    % Set up batch estimators
    batchBispEst = zeros(bispLen, batchNo);
    
    % Estimate batches
    if nargin==8
        for J=1:batchNo
            batchEst = zeros(bispLen, batchSampleSize(J));
            parfor M=batchLowerBounds(J):batchUpperBound(J)
                sampleSHC = image2shc(images(:, :, M), bandlimit, tDesign, interval, scalingParam);
                batchEst(:, M) = bispectrum(sampleSHC, bandlimit, CGs) ...
                    - sigma2*debiasingMatrix*cSHC2rSHC(sampleSHC);
            end
            batchBispEst(:, J) = mean(batchEst, 2);
        end
    elseif nargin==9
        for J=1:batchNo
            batchEst = zeros(bispLen, batchSampleSize(J));
            parfor M=batchLowerBounds(J):batchUpperBound(J)
                sampleSHC = image2shc(images(:, :, M), bandlimit, tDesign, interval, scalingParam, sh);
                batchEst(:, M) = bispectrum(sampleSHC, bandlimit, CGs) ...
                    - sigma2*debiasingMatrix*cSHC2rSHC(sampleSHC);
            end
            batchBispEst(:, J) = mean(batchEst, 2);
        end
    end
    
    % Estimate the bispectrum
    bispEst = mergeEstimatedBatches(batchBispEst, batchSampleSize);
else
    bispEst = zeros(bispLen, 1);
    if nargin==8
        for J=1:sampleSize
            sampleSHC = image2shc(images(:, :, J), bandlimit, tDesign, interval, scalingParam);
            sampleBisp = bispectrum(sampleSHC, bandlimit, CGs) ...
                    - sigma2*debiasingMatrix*cSHC2rSHC(sampleSHC);
            bispEst = (J-1)/J * bispEst + sampleBisp/J;
        end
    elseif nargin==9
        for J=1:sampleSize
            sampleSHC = image2shc(images(:, :, J), bandlimit, tDesign, interval, scalingParam, sh);
            sampleBisp = bispectrum(sampleSHC, bandlimit, CGs) ...
                    - sigma2*debiasingMatrix*cSHC2rSHC(sampleSHC);
            bispEst = (J-1)/J * bispEst + sampleBisp/J;
        end
    end
end