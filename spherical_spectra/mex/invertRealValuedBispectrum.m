function [invertedSHC, squaredResidual, output] ...
    = invertRealValuedBispectrum(bispectrum, bandlimit)
    
opts = optimoptions(@lsqnonlin, ...
    'SpecifyObjectiveGradient', true, ...
    'OptimalityTolerance', 10^-12, ...
    'FunctionTolerance', 10^-15, ...
    'StepTolerance', 10^-10, ...
    'Display', 'iter'); 
r = Inf;

x0 = rand((bandlimit+1)^2, 1);
%%
func = @(shc) inversionObjectiveFunc(shc, bispectrum, bandlimit);
% while r>10^-10
[invertedSHC, squaredResidual, ~, ~, output] = lsqnonlin(func, x0, [], [], opts);
r = sqrt(squaredResidual);
% end

end

%% The objective function
function [F, grad] = inversionObjectiveFunc(shc, bispectrum, L)
F = calculateBispectrumOfRealValuedFunction(shc, L) - bispectrum;
grad = calculateGradientOfBispectrumOfRealValuedFunction(shc, L)';
end