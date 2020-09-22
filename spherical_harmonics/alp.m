function out = alp(l, x)
% TODO: update docs, modernize the code and change the function name to
% something more readable.
% 
% Calculates the associated Legendre polynomial P_l^m at x for abs(m)<=l.
%   out(m+l+1, k) = P_l^m (x(k))

assert(real(l)==l && l>=0, 'l must be a non-negative integer.');

P1 = legendre(l, x); % Contains the ALP for m=0,1,2,...,l.

ms = (1:l)';
factors = (-1).^ms .* factorial(l-ms)./ factorial(l+ms);
P2 = factors .* P1(2:end, :);
P2 = flip(P2, 1);
out = [P2; P1];
