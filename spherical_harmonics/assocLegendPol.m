function alp = assocLegendPol(l, x)
%%
% Call format
%   alp = assocLegendPol(l, x)
% 
% Evaluate at x all the associated Legendre polynomial of order l. 
% 
% Input arguments
%   l       double      an integer, the order of the associated Legendre
%                       polynomial.
%   x       double      n-element array of numbers in [-1,1].
% 
% Output arguments
%   alp     double      (2l+1) x n array, such that
%                           out(m+l+1, k) = P_{l}^{m} (x(k))
%                       where P_{l}^{m} is the associated Legendre
%                       polynomial of order l and degree m (-l<=m<=l).
% 
% Notes
%   (1) The code performs no input checks.
%   (2) The code depends on MATLAB's alegendre function.
% 
% Reference
%   [1] Bateman, H. (1953). Higher Transcendental Functions (Vol. 1). 
%       McGraw-Hill Book Company.
% ***********************************************************
% Author    Ido Hadi
% Email     idohadi@mail.tau.ac.il
% Year      2020
% ***********************************************************

%% Compute associated Legendre polynomial
% Compute for non-negative degrees
alp = legendre(l, x);
% Compute for negative degress, using [1] (p. 140), eq. (7) and the fact
% Gamma(s+1)=s! for non-negative integer s
ms = (1:l)';
factors = (-1).^ms .* exp(gammaln(l-ms+1) - gammaln(l+ms+1));
alp = [flip(factors .* alp(2:end, :), 1); alp];
