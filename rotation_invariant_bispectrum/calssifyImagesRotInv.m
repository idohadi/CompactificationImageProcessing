function [avgedData, nearestNeighbors] = calssifyImagesRotInv(data, truncation, varargin)
%% TODO
% Call format
%   [avgedData, nearestNeighbors] = calssifyImagesRotInv(data, truncation, angularLimits)
%   [avgedData, nearestNeighbors] = calssifyImagesRotInv(data, truncation, angularLimits, __)
% 
% Identify close neighbors of images in data and denoise each image, using
% its neighbors.
% 
% 
% Input arguments
%   data                double      imageSize x imageSize x sampleSize
%                                   array, 
%   truncation          double      scalar, positive integer, bandlimit for
%                                   projection.
%   angularLimits       double 
% 
% Output arguments
%   avgedData           double      
%   nearestNeighbors    double      
% 
% Optional arguments
%   interval            Parameter for image2shc.
%   JaccardThreshold    The threshold over which we maintain an edge in the
%                       similarity matrix.
%   K                   Bispectrum denoising matrix.
%                       Default is false. In this case, the code computes
%                       it below.
%   Nneighbors          Number of nearest neighbors to find.
%   scalingParam        Projection scaling parameter to use in image2shc.
%   wpass               MATLAB's lowpass function parameter.
%                       If wpass=0, no low-pass filter is used.
% 
% Default optional arguments
%   JaccardThreshold    0.5
%   Nneighbors          50
%   wpass               0.05
% 
% Notes
%   None
% 
% Reference
%   [1] Halko, N., Martinsson, P.-G., Shkolnisky, Y., & Tygert, M. (2011). 
%       An Algorithm for the Principal Component Analysis of Large Data 
%       Sets. SIAM Journal on Scientific Computing, 33(5), 2580–2594. 
%       https://doi.org/10.1137/100804139
% ***********************************************************
% Author    Ido Hadi
% Email     idohadi@mail.tau.ac.il
% Year      2020
% ***********************************************************

%% Input handling
assert(size(data, 1)==size(data, 2), 'Images must be rectangular.');
imageSize = size(data, 1);
sampleSize = size(data, 3);
assert(isscalar(truncation) & truncation>=1 & round(truncation)==truncation, ...
    'Bandlimit must be a positive integer scalar.');

% Setting up optional input handling (name, value pairs)
p = inputParser;
addParameter(p, 'batchSize', 100, @(x) isscalar(x) & x>=1);
addParameter(p, 'JaccardThreshold', 0.2, @(x) isscalar(x) & x>=0);
addParameter(p, 'Nneighbors', 50, @(x) isscalar(x) & x>=1);
addParameter(p, 'wpass', 0.05, @(x) isscalar(x) & x>=0 & x<1);

% Process the optional input
parse(p, varargin{:});
batchSize = p.Results.batchSize;
JaccardThreshold = p.Results.JaccardThreshold;
Nneighbors = p.Results.Nneighbors;
wpass = p.Results.wpass;

k = 400;

%% Classify images
basis = ffb_basis(imageSize*ones(1, 2), truncation);
angularLimits = basis.k_max;

% Compute the length of the bispectrum vector
bLen = rotationInvariantBispectrumLength(truncation, angularLimits);

% Compute bispectra of data
clear bispectrum;
clear bispectrum_mex;

if wpass==0     % Don't use a low-pass filter
    rowFunc = struct('rowFunc', @rowFuncWODenoising, ...
            'args', {{data, basis, truncation, angularLimits, bLen}}, ...
            'm', sampleSize, ...
            'n', 2*bLen+angularLimits(1));
    [U, S, ~] = outOfCoreRandomizedSVD(rowFunc, k);
    b = U*S;
    clear U;
    clear S;
else            % Use a low-pass filter
    rowFunc = struct('rowFunc', @rowFuncWDenoising, ...
            'args', {{data, basis, truncation, angularLimits, bLen, wpass}}, ...
            'm', sampleSize, ...
            'n', 2*bLen+angularLimits(1));
    [U, S, ~] = outOfCoreRandomizedSVD(rowFunc, k);
    b = U*S;
    clear U;
    clear S;
end

% Compute nearest neighbors
[idx, D] = knnsearch(b, b, 'K', Nneighbors);

% Construct similarity matrix
W = sparse(repelem((1:sampleSize)', Nneighbors), ...
    reshape(idx', [numel(idx), 1]), 1, ...
    sampleSize, sampleSize, numel(idx));
W = W - diag(diag(W));

% Construct the Jaccard index to denoise the nearest neighbors graph
E = ones(size(W, 1), 1);
E = E*(E.');
JacInd = W*(W.')./( W*E + E*(W.') - W*(W.') );

% Denoise the similarity matrix
Wdenoised = W.*(W.'); % An edge a->b is kept only if there is an edge b->a
Wdenoised(JacInd<JaccardThreshold) = 0;

% Save the denoised similarity graph
nearestNeighbors = struct('W', W, 'D', D, 'Wdenoised', Wdenoised);

% Denoise data using the averaging
avgedData = zeros(size(data));
% TODO

end

%% Utility functions
function row = rowFuncWDenoising(N, data, basis, truncation, angularLimits, bLen, wpass)
% With denoising
coeff = basis.expand(lowpass(data(:, :, N), wpass));
row = rotationInvariantBispectrum(coeff, truncation, angularLimits, bLen);
end

function row = rowFuncWODenoising(N, data, basis, truncation, angularLimits, bLen)
% Without denoising
coeff = basis.expand(data(:, :, N));
row = rotationInvariantBispectrum(coeff, truncation, angularLimits, bLen);
end
