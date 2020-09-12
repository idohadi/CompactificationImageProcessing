function [invertedSHC, rootedResidual, output] = ibispectrum(b, bandlimit, tol)
%%
% Call format
%   invertedSHC = ibispectrum(b, bandlimit)
%   [invertedSHC, rootedResidual, output] = ibispectrum(b, bandlimit)
%   ibispectrum(b, bandlimit, tol)
% 
% Find a spherical harmonics coefficients vector approximately satisfying 
%   bispectrum(shc) = b.
% 
% Input arguments
%   b               double      bispectrum f spherical function
%                               represented by shc.
%   bandlimit       double      positive integer, the bandlimit of shc.
%   tol             double      tolerance for inversion loop
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
if nargin<=2
    tol = sqrt(10^-6);
end

%% Setup options
opts = optimoptions(@lsqnonlin, ...
    'SpecifyObjectiveGradient', true, ...
    'OptimalityTolerance', 10^-6, ...
    'FunctionTolerance', 10^-9, ...
    'StepTolerance', 10^-8, ...
    'MaxIterations', 300, ...
    'Display', 'iter'); 
rootedResidual = Inf;

%% Invert the bispectrum
nearOptimal = false;
nearOptimalTries = 0;
ef = 0;
while rootedResidual>tol
    if ~nearOptimal || nearOptimalTries>3
        nearOptimalTries = 0;
        initialSHC = randSHC(bandlimit);
%     initialSHC = r2c(initialSHC, bandlimit);
        initialSHC = cSHC2rSHC(initialSHC);
    else
        nearOptimalTries = nearOptimalTries + 1;
        direction = randn(size(invertedSHC));
        direction = direction/norm(direction);
        initialSHC = invertedSHC + 2*output.stepsize*direction;
        nearOptimal = false;
    end
    
    [invertedSHC, squaredResidual, ~, ef, output] = lsqnonlin(@inversionObjectiveFunc, initialSHC, [], [], opts);
    if output.firstorderopt<=10^-4 && squaredResidual <=10^-4
        nearOptimal = true;
    end
    rootedResidual = sqrt(squaredResidual);
    if ef==1 && rootedResidual<=5*10^-3
        break;
    end
end
invertedSHC = rSHC2cSHC(invertedSHC);

%% The objective function
function [F, grad] = inversionObjectiveFunc(shc)
[F, grad] = bispectrum(rSHC2cSHC(shc), bandlimit);
F = F - b;
end

end