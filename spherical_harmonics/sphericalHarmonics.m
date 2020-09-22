function sh = sphericalHarmonics(theta, phi, bandlimit)
%%
% Call format
%   sh = sphericalHarmonics(theta, phi, bandlimit)
% 
% Evaluate at (theta, phi) the spherical harmonics of order 0 to bandlimit 
% and all degress.
% 
% Convention
%   theta and phi represent points on the 2-sphere, as follows:
%   Variable    Interval
%   theta       [0,pi]          
%   phi         [0,2pi)         
% 
%   theta and phi satisfy
%   x(j) = cos(phi)*sin(theta)
%   y(j) = sin(phi)*sin(theta)
%   z(j) = cos(theta)
% 
%   Spherical harmonics of order l (0<=l<=bandlimit) and degree m (-l<=m<=l):
%   
%                                   ( 2l+1   (l-m)! )^(1/2)
%       Y_{l}^{m} (theta, phi) =    ( ---- * ------ )       * exp(i*m*phi) * P_{l}^{m} (cos(theta))
%                                   ( 4*pi   (l+m)! )
%   
%   where P_{l}^{m} (x) is the associated Legendre polynomial of order l
%   and degree m.
%       
%   TODO: make sure the code below actually computes this. There is a
%   chance I could save up on computing some of the factorials if I adopt
%   this (Varshalovic's) convention.
% 
% 
% Input arguments
%   theta       double      N-element array
%   phi         double      N-element array
%   bandlimit   double      positive integer
% 
% Output arguments
%   sh          double      N x (bandlimit+1)^2 array, such that
%                               shc(n, l*(l+1)+m+1) = Y_{l}^{m} (theta(n), phi(n))
% 
% Notes
%   (1) The code performs no input checks.
%   (2) TODO: refer to the book by Varshalovich about the convention of
%   spherical harmonics I used.
% 
% Reference
%   [1] Varshalovich, D. A., Moskalev, A. N., & Khersonskii, V. K. (1988). 
%       Quantum Theory of Angular Momentum. World Scientific. 
%       https://www.worldscientific.com/worldscibooks/10.1142/0270
% ***********************************************************
% Author    Ido Hadi
% Email     idohadi@mail.tau.ac.il
% Year      2020
% ***********************************************************

%% Compute spherical harmonics
theta = theta(:);
theta = theta.';
phi = phi(:);
phi = phi.';

fun = @(l, m) sqrt( (2*l+1)/(4*pi) * factorial(l-m)./factorial(l+m));

sh = cell(bandlimit, 1);
for l=1:bandlimit
    ms = (1:2*l+1) - l - 1;
    ms = ms(:);
    factors = fun(l, ms);
    
    sh{l} = ( exp(1i*ms.*phi) .* assocLegendPol(l, cos(theta)));
    sh{l} = ((-1).^ms) .* factors .* sh{l};
end
sh = vertcat(sh{:});
sh = [0.5/sqrt(pi) * ones(1, length(theta)); sh];
