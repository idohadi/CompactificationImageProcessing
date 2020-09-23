function [relativeDistance, alignedSHC2, optimalRotation, output] ...
        = alignSHC(shc1, shc2, bandlimit, tDesign, sequenceSize)
%%
% Call format
%   alignSHC(shc1, shc2, L, tDesign)
%   alignSHC(shc1, shc2, bandlimit, tDesign, sequenceSize)
%   relativeDistance = alignSHC(__)
%   [relativeDistance, alignedSHC2] = alignSHC(__)
%   [relativeDistance, alignedSHC2] = alignSHC(__)
%   [relativeDistance, alignedSHC2, optimalRotation] = alignSHC(__)
%   [relativeDistance, alignedSHC2, optimalRotation, output] = alignSHC(__)
%   
% Align two bandlimited spherical functions with bandlimit bandlimit,
% represented by spherical harmonics coefficients shc1 and shc2. The result
% is the spherical harmonics coefficients alignedSHC2 approximately 
% satisfying
%       min_{r is a 3d rotation} ||shc1 - rotatedSHC(shc2, r)||_2
% 
% The alignment relies on a brute force approach. The procedure follows the
% following stages:
% 1. Crude alignment. For every element r of a deterministic, 
%   low-discrepency sequence of rotations (see [1]) up to rotation with 
%   index sequenceSize, calculate the correlation of shc1 and 
%   rotateSHC(shc2, r). r1 is the rotation maximizing this correlation.
% 2. Refined alignment. Using MATLAB's fminsearch initialized at r1, find
%    a rotation maximizing said correlation.
% 
% Input arguments
%   shc1, shc2      double      (bandlimit+1)^2 x 1 arrays, complex vector 
%                               of spherical harmonics coefficients ordered 
%                               lexicographically.
%   bandlimit       double      positive integer, the bandlimit of shc1 and
%                               shc2.
%   tDesign         double      n x 3 array, a t-design. t is required to
%                               be greater than 2*bandlimit.
% 
% Optional arguments
%   Name                Default         Definition
%   sequenceSize        72*2^8          The maximum index in the sequence
%                                       of rotations to test in the crude
%                                       alignment stage.
% 
% Output arguments
%   relativeDistance double     scalar, norm(shc1 - alignedSHC2)/norm(shc1).
%   alignedSHC2      double     (bandlimit+1)^2 x 1 array, shc2
%                               approximately aligned to shc1 in the sense
%                               discussed above.
%   optimalRotation  double     4 x 1 array, the optimal rotation in
%                               quaternion representation. That is, 
%                                   alignedSHC2  = rotateSHC(shc2, optimalRotation)
%   output           double     struct, containing various output arguments
%                               from MATLAB's fminsearch and the two stages.
%                               The fminsearch output is documented in 
%                               MATLAB's docs.
% 
% Notes
%   (1) This function performs no input checks.
%   (2) This function depends on the SmallRotationToolbox repository. This
%       repository contains an implementation of the algorithm of [1],
%       which is necessary for the crude alignment stage.
%       https://github.com/idohadi/SmallRotationToolbox
% 
% Reference
%   [1] Yershova, A., Jain, S., LaValle, S. M., & Mitchell, J. C. (2009). 
%       Generating Uniform Incremental Grids on SO(3) Using the Hopf 
%       Fibration. The International Journal of Robotics Research, 29(7), 
%       801–812. https://doi.org/10.1177/0278364909352700
% ***********************************************************
% Author    Ido Hadi
% Email     idohadi@mail.tau.ac.il
% Year      2020
% ***********************************************************

%% Input validation
narginchk(4, 5);

if nargin == 4
    sequenceSize = 72*2^8;
end

%% Align SHCs
output = struct();

% Stage 1: crude alignment
corrs = zeros(sequenceSize, 1);
parfor r=1:sequenceSize
    rotation = quasiRandRotation(r);
    rotatedSHC2 = rotateSHC(shc2, bandlimit, rotation, tDesign);
    corrs(r) = corr([real(shc1); imag(shc1)], ...
        [real(rotatedSHC2); imag(rotatedSHC2)]);
end
[crudeMaxCorrelation, maxInd] = max(corrs);

output.crudeMaxCorrelation = crudeMaxCorrelation;
output.crudeMaxCorrelationIndex = maxInd;

% Stage 2: refined alignment
[theta, phi, ~] = cart2sph2(tDesign(:, 1), tDesign(:, 2), tDesign(:, 3));
sh = sphericalHarmonics(theta, phi, bandlimit);
shc1tDesignVec = sh.' * shc1;
shc1tDesignVec = [real(shc1tDesignVec); imag(shc1tDesignVec)];

initPoint = quasiRandRotation(maxInd);
opts = optimset('Display', 'off', ...
                    'TolFun', 10^-10, ...
                    'TolX', 10^-10, ...
                    'MaxFunEvals', 400);
[optimalRotation, refinedMaxCorrelation, exitflag, fminsearchOutput] ...
            = fminsearch(@(x) objFunc(x), initPoint, opts);

optimalRotation = optimalRotation/norm(optimalRotation);
if optimalRotation(1)<0
    optimalRotation = - optimalRotation;
end

alignedSHC2 = rotateSHC(shc2, bandlimit, optimalRotation, tDesign);
relativeDistance = norm(shc1 - alignedSHC2, 2)/norm(shc1, 2);

output.refinedMaxCorrelation = 1 - refinedMaxCorrelation;
output.exitflag = exitflag;
output.fminsearchOutput = fminsearchOutput;

% fminsearch objective function
function out = objFunc(x)
    x = x/norm(x);
    if x(1)<0
        x = -x;
    end

    shc2tDesignVec = sh.' * rotateSHC(shc2, bandlimit, x, tDesign);
    shc2tDesignVec = [real(shc2tDesignVec); imag(shc2tDesignVec)];

    out = corr(shc1tDesignVec, shc2tDesignVec);
    out = 1 - out;
end

end
