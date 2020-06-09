function grad = bispgrad(shc, L)
grad = calculateGradientOfBispectrumOfRealValuedFunction(shc, L);
grad(isnan(grad)) = 0;