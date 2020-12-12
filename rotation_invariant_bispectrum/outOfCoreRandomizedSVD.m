function [U, S, V] = outOfCoreRandomizedSVD(rowFunc, m, n, k, varagin)
%% TODO
% Call format
%   [U, S, V] = outOfCoreRandomizedSVD(rowFunc, k)
%   [U, S, V] = outOfCoreRandomizedSVD(rowFunc, k, __)
% 
% Identify close neighbors of images in data and denoise each image, using
% its neighbors.
% 
% 
% TODO edit docs below
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
%   (1) Notations follow [1].
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
% TODO
l = k + 2;
i = 2;
batchSize = 100;

func = @(x) rowFunc.fucn(x, rowFunc.args{:});

%% Compute the randomized SVD
% Generate batch limits
batches = [0, batchSize:batchSize:m];
if batches(end)<m
    batches(end+1) = m;
end

% Alg: step 1
G = randn(n, l);

% Initialize H
H = zeros(m, (i+1)*l);

% Compute H^0
H(:, 1:l) = AprodB(G);

% Compute H^s for  s=1, ..., i
for I=1:i
    H(:, l*(I-1)+1:l*I) = AprodB(AtprodB(H(:, l*(I-2)+1:l*(I-1))));
end

% Alg: step 2
[Q, ~, ~] = qr(H);

% Alg: step 3
T = AtprodB(Q);

% Alg: step 4
[Vtilde, Stilde, W] = svd(T, 'econ');

% Alg: step 5
Utilde = Q*W;

% Step 6
U = Utilde(:, 1:k);
V = Vtilde(:, 1:k);
S = Stilde(1:k, 1:k);


%% Utility functions
function P = AprodB(B)
% Compute the product P = A * B
assert(size(B, 1)==n, ...
    'Number of rows of B must equal number of columns of A.');

P = zeros(m, size(B, 2));

for b=1:length(batches)-1
    A = zeros(batches(b+1) - batches(b), n);
    batchLowLim = batches(b);
    parfor J=batches(b)+1:batches(b+1)
        A(J - batchLowLim, :) = func(J);
    end
    P(batches(b)+1:batches(b+1), :) = A*B;
end
end

function P = AtprodB(B)
% Compute the product P = A.' * B
assert(size(B, 1)==m, ...
    'Number of rows of B must equal number of rows of A.');

P = zeros(n, size(B, 2));

for b=1:length(batches)-1
    At = zeros(n, batches(b+1) - batches(b));
    batchLowLim = batches(b);
    parfor J=batches(b)+1:batches(b+1)
        At(:, J - batchLowLim) = func(J);
    end
    P = P + At .* B(batches(b)+1:batches(b+1), :);
end
end

end