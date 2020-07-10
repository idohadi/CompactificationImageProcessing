function F = inversionObjectiveFuncManopt(shc, bispectrum, L, bispInds, symComp)
    shc = struct2vecManopt(shc, L);
    shc = r2c(shc, L);
    shc = SHCProdFactor(shc, 1/sqrt(2), L);
    
    [F, ~] = inversionObjectiveFuncReal(shc, bispectrum, L, bispInds, symComp);
    F = norm(F, 2)^2;
end
