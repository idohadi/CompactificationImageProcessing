function b = calculateBispectrum_M(shc, L)
b = calculateBispectrum_MATMIM(shc, conj(shc), L);
