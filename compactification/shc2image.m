function [im, sh] = shc2image(shc, bandlimit, imageSize, interval, scalingParam, sh)
%%
% Call format
%   im = shc2image(shc, bandlimit, tDesign, interval, scalingParam)
%   im = shc2image(shc, bandlimit, tDesign, interval, scalingParam, sh)
%   [im, sh] = shc2image(__)
% 
% Creates an image by projecting a bandlimited spherical function into R^2.
% 
% The image is assumed to be the values of a function on an equidistant 
% grid in the cube interval^2.
% 
% Input arguments
%   shc             double      (bandlimit+1)^2 x 1 array, complex vector 
%                               of spherical harmonics coefficients ordered 
%                               lexicographically.
%   bandlimit       double      positive integer, the bandlimit of shc.
%   imageSize       double      im will be of size imageSize x imageSize.
%   interval        double      1 x 2 array, interval(1) is the lower bound
%                               of the interval and interval(2).
%   scalingParam    double      positive number, scaling parameter for 
%                               the projection.
% 
% Output arguments
%   im              double      n x n array, an image
% 
% Notes
%   This function performs no input checks.
% 
% Reference
%   None
% ***********************************************************
% Author    Ido Hadi
% Email     idohadi@mail.tau.ac.il
% Year      2020
% ***********************************************************

%% Compute the image
% Generate a grid that the image is evaluated at
X = linspace(interval(1), interval(2), imageSize);
[X, Y] = meshgrid(X, X);
X = X(:);
Y = Y(:);

% Project grid onto sphere
[R2phi, R2rho] = cart2pol2(X, Y);
[theta, phi] = KondorProj(R2phi, R2rho, scalingParam);

% Evaluate image
if nargin<6
sh = sphericalHarmonics(theta, phi, bandlimit);
end
im = sh.'*shc;
im = reshape(im, [imageSize, imageSize]);
