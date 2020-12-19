function [U, S, V] = outOfCoreRandomizedSVD(A, k, varargin)
%% TODO
% Call format
%   [U, S, V] = outOfCoreRandomizedSVD(A, k)
%   [U, S, V] = outOfCoreRandomizedSVD(A, k, __)
%   [U, S, V] = outOfCoreRandomizedSVD(rowFunc, k)
%   [U, S, V] = outOfCoreRandomizedSVD(rowFunc, k __)
% 
% Compute a rank k approximation of a matrix using a randomized SVD 
% algorithm. U*S*V' is a k rank matrix appproximating the input matrix.
% 
% A can either be an m x n array or a structure. If it is a structure, it must
% have the following fields:
%   m       -   the number of rows of the matrix.
%   n       -   the number of columns of the matrix.
%   rowFunc -   a function handle such that rowFunc(J, args{:}) returns the
%               J-th row of the matrix.
%   args    -   auxilary arguments for rowFunc.
% 
% 
% Input arguments
%   A       double      an input matrix, specified in one of the two
%                       formats presented above.
%   k       double      the output is to be a rank k approximation of A.
% 
% Output arguments
%   U, S, V double      matrices such that 
%                           U*S*V' is approximately a rank k approximation
%                       of the input matrix A.
% 
% Optional arguments
%   l, i        Utility parameters. See [1].
%   batchSize   Batch sizes for parallel computation. Relevant only when
%               the input matrix is given as a row-generating function.
% 
% Default optional arguments
%   l           k + 2
%   i           2
%   batchSize   100
% 
% Notes
%   (1) This is an implementation of the algorithm presented in [1]. This
%       implementation assumes the input matrix is given in one of two
%       forms: either the matrix itself is given, or a row-generating
%       function is given.
%   (2)  Notations follow [1].
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
p = inputParser;

if ismatrix(A) && isnumeric(A)
    m = size(A, 1);
    n = size(A, 2);
    
    structFlag = false;
elseif isstruct(A)
    assert(isfield(A, 'm'), 'm (number of rows) must be a field.');
    assert(isfield(A, 'n'), 'n (number of columns) must be a field.');
    assert(isfield(A, 'rowFunc'), 'rowFunc (row-generating function) must be a field.');
    assert(isfield(A, 'args'), 'args (row-generating function auxilary arguments) must be a field.');
    
    m = A.m;
    n = A.n;
    func = @(x) A.rowFunc(x, A.args{:});
    
    structFlag = true;
else
    error('Unidentified input matrix format.');
end

% Process the optional input
addParameter(p, 'batchSize', 100, @(x) isscalar(x) & x>=1);
addParameter(p, 'l', k+2, @(x) isscalar(x) & x>=k);
addParameter(p, 'i', 2, @(x) isscalar(x) & x>=0);

parse(p, varargin{:});
batchSize = p.Results.batchSize;
i = p.Results.i;
l = p.Results.l;


%% Compute the randomized SVD
if structFlag
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
        H(:, l*(I-1)+1:l*I) = AprodB(AtprodB(H(:, l*(I-1)+1:l*I)));
    end

    % Alg: step 2
    [Q, ~, ~] = qr(H, 0);

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
else
    % Alg: step 1
    G = randn(n, l);

    % Initialize H
    H = zeros(m, (i+1)*l);

    % Compute H^0
    H(:, 1:l) = A*G;

    % Compute H^s for  s=1, ..., i
    for I=1:i
        H(:, l*(I-1)+1:l*I) = A*(A.'*H(:, l*(I-1)+1:l*I));
    end

    % Alg: step 2
    [Q, ~, ~] = qr(H, 0);

    % Alg: step 3
    T = A.'*Q;

    % Alg: step 4
    [Vtilde, Stilde, W] = svd(T, 'econ');

    % Alg: step 5
    Utilde = Q*W;

    % Step 6
    U = Utilde(:, 1:k);
    V = Vtilde(:, 1:k);
    S = Stilde(1:k, 1:k);
end

%% Utility functions
function P = AprodB(B)
% Compute the product P = A * B
assert(size(B, 1)==n, ...
    'Number of rows of B must equal number of columns of A.');

P = zeros(m, size(B, 2));

for b=1:length(batches)-1
    Arows = zeros(batches(b+1) - batches(b), n);
    batchLowLim = batches(b);
    parfor J=batches(b)+1:batches(b+1)
        Arows(J - batchLowLim, :) = feval(func, J);
    end
    P(batches(b)+1:batches(b+1), :) = Arows*B;
end
end

function P = AtprodB(B)
% Compute the product P = A.' * B
assert(size(B, 1)==m, ...
    'Number of rows of B must equal number of rows of A.');

P = zeros(n, size(B, 2));

for b=1:length(batches)-1
    Arowst = zeros(n, batches(b+1) - batches(b));
    batchLowLim = batches(b);
    parfor J=batches(b)+1:batches(b+1)
        Arowst(:, J - batchLowLim) = feval(func, J);
    end
    P = P + Arowst * B(batches(b)+1:batches(b+1), :);
end
end

end