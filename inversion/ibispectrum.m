function [invertedSHC, rootedResidual, output] = ibispectrum(b, bandlimit, x0)
%%
% Call format
%   ibispectrum(b, bandlimit)
%   ibispectrum(b, bandlimit, x0)
%   invertedSHC = ibispectrum(__)
%   [invertedSHC, rootedResidual] = ibispectrum(__)
%   [invertedSHC, rootedResidual, output] = ibispectrum(__)
% 
% Find a spherical harmonics coefficients vector approximately satisfying 
%   bispectrum(shc) = b.
% 
% Input arguments
%   b               double      bispectrum f spherical function
%                               represented by shc.
%   bandlimit       double      positive integer, the bandlimit of shc.
%   x0              double      (bandlimit+1)^2 complex array, an inital
%                               guess.
% 
% Output arguments
%   invertedSHC     double      complex vector array approximately 
%                               satisfying bispectrum(invertedSHC) = b.
%   rootedResidual  double      The 2-norm distance between b and
%                               bispectrum(invertedSHC). That is,
%                                   ||b - bipsectrum(invertedSHC)||_2 .
%   output          struct      Various parameters about the optimization
%                               process performed by MATLAB's lsqnonlin.
%                               for details, see docs of lsqnonlin.
% 
% Notes
%   (1) This function depends on bispectrum_mex.
%   (2) This function inverts the given bispectrum by iteratively solving
%       an optimization problem.
% 
% Reference
%   None
% ***********************************************************
% Author    Ido Hadi
% Email     idohadi@mail.tau.ac.il
% Year      2020
% ***********************************************************

%% Input validation
assert(bandlimit==round(bandlimit) && bandlimit>0, ...
    'bandlimit must be a positive integer.');

%% Setup options
pow = -8;
opts = optimoptions(@lsqnonlin, ...
    'SpecifyObjectiveGradient', true, ...
    'OptimalityTolerance', 10^pow, ...
    'FunctionTolerance', 10^pow, ...
    'StepTolerance', 10^pow, ...
    'Display', 'off'); 

%% Invert the bispectrum
if nargin<3
    initialSHC = cSHC2rSHC(randSHC(bandlimit));
else
    initialSHC = cSHC2rSHC(x0);
end

[invertedSHC, squaredResidual, ~, ~, output] ...
    = lsqnonlin(@inversionObjectiveFunc, initialSHC, [], [], opts);
rootedResidual = sqrt(squaredResidual);
invertedSHC = rSHC2cSHC(invertedSHC);

%% The objective function
function [F, grad] = inversionObjectiveFunc(shc)
[F, grad] = bispectrum(rSHC2cSHC(shc), bandlimit);
F = F - b;
end

end