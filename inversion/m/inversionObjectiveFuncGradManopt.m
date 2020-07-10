function grad = inversionObjectiveFuncGradManopt(shc, bispectrum, L, bispInds, symComp)
    shc = struct2vecManopt(shc, L);
    shc = r2c(shc, L);
    shc = SHCProdFactor(shc, 1/sqrt(2), L);
    
    [F, grad] = inversionObjectiveFuncReal(shc, bispectrum, L, bispInds, symComp);
    
    grad = 2*grad'*F;
    grad = vec2structManopt(grad, L);
end
