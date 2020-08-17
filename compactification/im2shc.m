function shc = im2shc(im, evalPoints, N, assocLegMat, omega, phi)
%%
% Call format
%   shc = im2shc(im, evalPoints, N)
% 
% Interpolates the image im to certain evaluation points and use them to
% estimate the spherical harmonics coefficients of the projection of the
% image onto the sphere.
% 
% Input arguments
%   im              double      n x m array, an image.
%   evalPoints      double      k x 2 array, points in R^2.
%   N               double      positive integer, the bandlimit of shc is
%                               2*N-1.
%   assocLegMat     double      cell array of matrices such that
%                               assocLegMat{m+2*N} is a (2*M-1) x (4*N-m)
%                               array satisfying
%                                   assocLegMat{m+2*N}(i, l-m+1) 
%                                       = >P_{l}^{|m|} (cos(theta(i)))
%                                   for 
%                                       0<=m<=2N-1 and m<=l<=2*N-1,
%                               where >P_{l}^{|m|} is a scaled, normalized
%                               associated Legendre function.
%   omega           double      2N x 1 array, Guass-Legendre cubature
%                               weights.
%   theta           double      2N x 1 array, arccos of Guass-Legendre 
%                               cubature points.
%   phi             double      (4N-1) x 1 array, grid on [0, 2pi].
%   
% Output arguments
%   shc             double      (2N)^2 x 1 array, complex vector of
%                               spherical harmonics coefficients ordered 
%                               lexicographically.
% 
% Notes
%   This function depends on the repository
%   FastSphericalHarmonicsTransform.
% 
% Reference
%   None
% ***********************************************************
% Author    Ido Hadi
% Email     idohadi@mail.tau.ac.il
% Year      2020
% ***********************************************************
%% Compute the sphericah harmonics coefficients
% Compute grid on which such that im(k) is the value given in (X(k), Y(k))
X = linspace(-0.5, 0.5, size(im, 2));
Y = linspace(-0.5, 0.5, size(im, 1));
[X, Y] = meshgrid(X, Y);

% Interpolate the given values to evalPoints
interpVals = interp2(X, Y, im, ...
                    evalPoints(:, 1), evalPoints(:, 2), ...
                    'cubic', 0);
interpVals = reshape(interpVals, [2*N, 4*N-1]);

% Estimate the spherical harmonics coefficients
shc = isht(interpVals, N, assocLegMat, omega, phi);
