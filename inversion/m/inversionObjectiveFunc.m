function [F, grad] = inversionObjectiveFunc(shc, bispectrum, L)
    if nargout==1
        F = calculateBispectrum_MATMIM(shc, conj(shc), L) - bispectrum;
    else
        [F, grad] = calculateBispectrum_MATMIM(shc, conj(shc), L);
        F = F - bispectrum;
        grad = grad.';
    end
end
