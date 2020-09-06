function [invertedSHC, rootedResidual, output] ...
    = invertComplexValuedBispectrum(bispectrum, bandlimit, varargin)
%TODO: Docs

%% Define parameters
opts = optimoptions(@lsqnonlin, ...
    'SpecifyObjectiveGradient', true, ...
    'OptimalityTolerance', 10^-12, ...
    'FunctionTolerance', 10^-15, ...
    'StepTolerance', 10^-10, ...
    'Display', 'iter', ...
    'CheckGradients', true); 
r = Inf;

tol = 10^-8;
if nargin>2
    tol = varargin{1};
end
    
%% Inversion loop
func = @(shc) inversionObjectiveFunc(cfy(shc), bispectrum, bandlimit);
while r>tol
%     x0 = randomNormalizedSHC(bandlimit, 1);
    x0 = randomNormalizedSHC(bandlimit, 0);
    x0 = rlfy(r2c(x0, bandlimit));
    [invertedSHC, rootedResidual, ~, ~, output] = lsqnonlin(func, x0, [], [], opts);
    r = output.firstorderopt;
end
invertedSHC = cfy(invertedSHC);
rootedResidual = sqrt(rootedResidual);

end

%% The objective function
function [F, grad] = inversionObjectiveFunc(shc, bispectrum, L)
    if nargout==1
        F = calculateBispectrum_MATMIM(shc, conj(shc), L) - bispectrum;
    else
        [F, grad] = calculateBispectrum_MATMIM(shc, conj(shc), L);
        F = F - bispectrum;
        grad = grad.';
    end
end