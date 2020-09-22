function out = evalYt(theta, phi, t)
% TODO: update docs, modernize the code and change the function name to
% something more readable.
% 
% out = shY(theta, phi, t)     evaluates the matrix Y_t (X_N)
%                               for t-design, as defined in
%                                   AN el al, 2010, eq. (2.3).
%                               X_N is reperesnted using 
%                               spherical coordinates.
%                               theta, phi must be vectors of 
%                               identical sizes.
% 
% NOTE: 
%   Spherical coordinates were defined in
%       Varshalovich et al, p. 4, eq. (8).
%   Spherical harmonics were defined in 
%       Varshalovich et al, p. 133, eq. (1). 
%   Y^0_t matrix were defined in
%       An et al, 2010, p. 2137, eq. (2.3).
% 
% REFERENCES:
%   AN, C., Chen, X., Sloan, I., & Womersley, R. (2010). 
%       Well Conditioned Spherical Designs for Integration 
%       and Interpolation on the Two-Sphere. 
%       SIAM Journal on Numerical Analysis, 48(6), 2135�2157. 
%       https://doi.org/10.1137/100795140
%   Varshalovich, D. A., Moskalev, A. N., & Khersonskii, V. K. (1988). 
%       Quantum Theory of Angular Momentum. World Scientific. 
%       Retrieved from https://www.worldscientific.com/worldscibooks/10.1142/0270
% 

narginchk(3, 3);
% assert(all(theta>=0) && all(theta<=pi), 'theta must be between 0 and pi.');
% assert(all(phi>=0) && all(phi<=2*pi), 'phi must be between 0 and 2pi.');
% assert(t>=1 && round(t)==t, 't must be a positive integer.');

theta = theta(:);
theta = theta.';
phi = phi(:);
phi = phi.';

fun = @(l, m) sqrt( (2*l+1)/(4*pi) * factorial(l-m)./factorial(l+m));

out = cell(t, 1);
for l=1:t
    ms = (1:2*l+1) - l - 1;
    ms = ms(:);
    factors = fun(l, ms);
    
    out{l} = ( exp(1i*ms.*phi) .* alp(l, cos(theta)));
    out{l} = ((-1).^ms) .* factors .* out{l};
end
out = vertcat(out{:});
out = [0.5/sqrt(pi) * ones(1, length(theta)); out];
