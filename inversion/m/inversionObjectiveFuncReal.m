function [F, grad] = inversionObjectiveFuncReal(shc, bispectrum, L, bispInds, symComp)
    [F, grad] = inversionObjectiveFunc(shc, bispectrum, L);
    F = F(bispInds);
    grad = grad(bispInds, :)*symComp;
end
