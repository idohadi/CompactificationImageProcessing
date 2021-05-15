function [G, D] = classifyImages(data, bandlimit, varargin)
%%
% Call format
%   [G, D] = classifyImages(data, bandlimit)
%   [G, D] = classifyImages(data, bandlimit, __)
% 
% Identify close neighbors of images in data and denoise each image, using
% its neighbors.
% 
% 
% Input arguments
%   data                double      imageSize x imageSize x sampleSize
%                                   array, 
%   bandlimit           double      scalar, positive integer, bandlimit for
%                                   projection.
% 
% Output arguments
%   G                   double      the nearest neighbors graph in a
%                                   matrix.
%   D                   double      The distances from MATLAB's nearest
%                                   neighbors function.
% 
% Optional arguments
%   interval            Parameter for image2shc.
%   K                   Bispectrum denoising matrix.
%                       Default is false. In this case, the code computes
%                       it below.
%   lowRank             false if a low-rank approximation is not to be used.
%                       Otherwise, the rank of the lonk rank to be used.
%                       otherwise.
%   Nneighbors          Number of nearest neighbors to find.
%   scalingParam        Projection scaling parameter to use in image2shc.
%   wpass               MATLAB's lowpass function parameter.
%                       If wpass=0, no low-pass filter is used.
% 
% Default optional arguments
%   interval            [-0.5, 0.5]
%   K                   false
%   lowRank             false
%   Nneighbors          50
%   scalingParam        1.5
%   wpass               0.05
% 
% Notes
%   None
% 
% Reference
%   None

% ***********************************************************
% Author    Ido Hadi
% Email     idohadi@mail.tau.ac.il
% Year      2020
% ***********************************************************

%% Input handling
assert(size(data, 1)==size(data, 2), 'Images must be rectangular.');
imageSize = size(data, 1);
sampleSize = size(data, 3);
assert(isscalar(bandlimit) & bandlimit>=1 & round(bandlimit)==bandlimit, ...
    'Bandlimit must be a positive integer scalar.');

% Setting up optional input handling (name, value pairs)
p = inputParser;
addParameter(p, 'interval', [-0.5, 0.5], @(x) numel(x)==2 & x(1)<x(2));
addParameter(p, 'K', false);
addParameter(p, 'lowRank', false, @(x) x==false | (isscalar(x) & x>=1 & round(x)==x));
addParameter(p, 'Nneighbors', 50, @(x) isscalar(x) & x>=1);
addParameter(p, 'scalingParam', 1.5, @(x) isscalar(x) & x>0);
addParameter(p, 'sigma2', 1, @(x) isscalar(x) & x>0);
addParameter(p, 'wpass', 0.05, @(x) isscalar(x) & x>=0 & x<1);


% Process the optional input
parse(p, varargin{:});
interval = p.Results.interval;
K = p.Results.K;
lowRank = p.Results.lowRank;
Nneighbors = p.Results.Nneighbors;
scalingParam = p.Results.scalingParam;
sigma2 = p.Results.sigma2;
wpass = p.Results.wpass;

%% Classify images
% Load t-design
td = loadtd(2*bandlimit);

% Declare global CGs variable
global CGs;

% Generating debiasing matrices
if K==false
    U = buildU(bandlimit, imageSize, td, interval, scalingParam);
    K = buildK(bandlimit, U);
end

% Compute the length of the bispectrum vector
shc = image2shc(data(:, :, 1), bandlimit, td, interval, scalingParam);
b = bispectrum(shc, bandlimit, CGs);
bLen = size(b, 1);
clear shc;
clear b;

% Compute bispectra of data
b = zeros(sampleSize, bLen);
clear bispectrum;
clear bispectrum_mex;

if wpass==0
    % Don't use a low-pass filter
    parfor n=1:sampleSize
        shc = image2shc(data(:, :, n), bandlimit, td, interval, scalingParam);
        b(n, :) = bispectrum(shc, bandlimit, CGs) - sigma2*K*cSHC2rSHC(shc);
    end
else
    % Use a low-pass filter
    parfor n=1:sampleSize
        shc = image2shc(lowpass(data(:, :, n), wpass), bandlimit, td, interval, scalingParam);
        b(n, :) = bispectrum(shc, bandlimit, CGs) - sigma2*K*cSHC2rSHC(shc);
    end
end

% Compute rank approximation of b
if ~islogical(lowRank)
    [U, S, ~] = outOfCoreRandomizedSVD(b, lowRank);
    b = U*S;
end

% Compute nearest neighbors
[idx, D] = knnsearch(b, b, 'K', Nneighbors);

% Construct similarity matrix
G = sparse(repelem((1:sampleSize)', Nneighbors), ...
    reshape(idx', [numel(idx), 1]), 1, ...
    sampleSize, sampleSize, numel(idx));
G = spdiags(zeros(sampleSize, 1), 0, G);
