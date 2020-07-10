function [F, grad] = inversionObjectiveFuncRealFM(shc, bispectrum, L, K, bispInds, symComp, bispIndsFM)
% shc needs to be a K-bandlimited (K+1)^2 length array of SHCs of a 
% real-valued function

shc = r2c(shc, K);
shc = [shc; 2*(rand((L+1)^2-(K+1)^2, 1) - 1) + 1i*2*(rand((L+1)^2-(K+1)^2, 1) - 1)];
[F, grad] = inversionObjectiveFunc(shc, bispectrum, L);
F = F(bispInds & bispIndsFM(K+1, :));
grad = grad(bispInds & bispIndsFM(K+1, :), :)*symComp;
grad = grad(:, (K^2+1):(K+1)^2);
