% Align spherical harmonics. the refinement step should be done here. The crude step in mex function.
% TODO: write docs and function itself
% 

% The following is from an old codebase
function [relativeDistance, alignedSHC2, optimalRotation, ...
    refinedMaxCorrelation, crudeMaxCorrelation, fminserachOutput] ...
        = alignSHC(shc1, shc2, L, tDesign, varargin)
%% Handle optional arguments
p = inputParser;
addParameter(p, 'SequenceSize', 2^8*72);
parse(p, varargin{:});
sequenceSize = p.Results.SequenceSize;

%% Evaluate Yt on uniform points on sphere
[theta, phi] = sphcart2sph(tDesign(:, 1), tDesign(:, 2), tDesign(:, 3));
Yt = evalYt(theta, phi, L);

%% First crude optimization
corrs = zeros(sequenceSize, 1);
parfor r=0:(sequenceSize-1)
    rotation = generateUniformRotationSequencePoint(r);
    rotatedSHC2 = rotateSphericalHarmonicsByEstimation(shc2, L, tDesign, rotation);
    corrs(r+1) = corr([real(shc1); imag(shc1)], ...
        [real(rotatedSHC2); imag(rotatedSHC2)]);
end
[crudeMaxCorrelation, maxInd] = max(corrs);

%% Refining the optimization using optimziation function
initPoint = generateUniformRotationSequencePoint(maxInd-1);
options = optimset('Display', 'off', ...
                    'TolFun', 10^-10, ...
                    'TolX', 10^-10, ...
                    'MaxFunEvals', 400);
[optimalRotation, refinedMaxCorrelation, ~, fminserachOutput] ...
            = fminsearch(@(x) objFunc(x, shc1, shc2, Yt, L, tDesign), ...
                initPoint, options);

optimalRotation = optimalRotation/norm(optimalRotation);
if optimalRotation(2)<0
    optimalRotation = - optimalRotation;
end

% optimalRotation = generateUniformRotationSequencePoint(maxInd-1);
alignedSHC2 = rotateSphericalHarmonicsByEstimation(shc2, L, tDesign, optimalRotation);
relativeDistance = norm(shc1 - alignedSHC2, 2)/norm(shc1, 2);

end

%% Objective function
function out = objFunc(x, shc1, shc2, Yt, L, tDesign)
x = x/norm(x);
if x(2)<0
    x = -x;
end

out = corr(real(Yt.' * shc1), ...
            real(Yt.' * rotateSphericalHarmonicsByEstimation(shc2, L, tDesign, x)));
out = 1 - out;
end