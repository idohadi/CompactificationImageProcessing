function [invertedSHC, rootedResidual, output] ...
    = invertComplexValuedBispectrum(bispectrum, bandlimit)
%TODO: Docs

%% Define parameters
opts = optimoptions(@lsqnonlin, ...
    'SpecifyObjectiveGradient', true, ...
    'OptimalityTolerance', 10^-12, ...
    'FunctionTolerance', 10^-15, ...
    'StepTolerance', 10^-10, ...
    'Display', 'iter', ...
    'CheckGradients', true); 

%% Inversion loop
func = @(shc) inversionObjectiveFunc(cfy(shc), bispectrum, bandlimit);
% while r>10^-10
x0 = randomNormalizedSHC(bandlimit, 1);
[invertedSHC, rootedResidual, ~, ~, output] = lsqnonlin(func, x0, [], [], opts);
invertedSHC = cfy(invertedSHC);
rootedResidual = sqrt(rootedResidual);
% end

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