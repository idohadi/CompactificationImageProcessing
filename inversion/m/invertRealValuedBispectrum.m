function [invertedSHC, rootedResidual, output] = invertRealValuedBispectrum(bispectrum, bandlimit)
%%
% opts = optimoptions(@lsqnonlin, ...
%    'SpecifyObjectiveGradient', true, ...
opts = optimoptions(@lsqnonlin, ...
    'OptimalityTolerance', 10^-12, ...
    'FunctionTolerance', 10^-15, ...
    'StepTolerance', 10^-10, ...
    'Display', 'iter'); 


%%
func = @(shc) inversionObjectiveFunc(shc, bispectrum, bandlimit);

c = 0;
r = Inf;
% while (c<=100 || r<=10^-8)
    x0 = randomNormalizedSHC(bandlimit, 0);
    [invertedSHC, rootedResidual, ~, ~, output] = lsqnonlin(func, x0, [], [], opts);
    r = output.firstorderopt;
    c = c + 1;
% end
rootedResidual = sqrt(rootedResidual);

end

%% The objective function
function F = inversionObjectiveFunc(shc, bispectrum, L)
[F, grad] = calculateBispectrum(shc, L, 0);
F = F - bispectrum;
grad = grad.';
end