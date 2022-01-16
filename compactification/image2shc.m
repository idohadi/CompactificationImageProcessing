function [shc, sh] = image2shc(im, bandlimit, tDesign, interval, scalingParam, sh)
%%
% Call format
%   shc = image2shc(im, bandlimit, tDesign, interval, scalingParam)
%   shc = image2shc(im, bandlimit, tDesign, interval, scalingParam, sh)
%   [shc, sh] = image2shc(__)
%   
% Interpolates the image im to a select t-design points and use them to
% estimate the spherical harmonics coefficients of the projection of the
% image onto the sphere. 
% 
% The image is assumed to be the values of a function on an equidistant 
% grid in the cube interval^2.
% 
% Input arguments
%   im              double      n x n array, an image
%   bandlimit       double      positive integer, the bandlimit of shc
%   tDesign         double      N x 3 array, a spherical design (a 
%                               t-design) in Cartesian coordiantes.
%   interval        double      1 x 2 array, interval(1) is the lower bound
%                               of the interval and interval(2).
%   scalingParam    double      positive number, scaling parameter for 
%                               the projection.
% 
% Output arguments
%   shc             double      (bandlimit+1)^2 x 1 array, complex vector 
%                               of spherical harmonics coefficients ordered 
%                               lexicographically.
% 
% Optional input/output arguments
%   sh              double      Spherical harmonics matrix. If it is given
%                               as an input, it is not recomputed below.
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

%% Compute the sphericah harmonics coefficients
% t-design in spherical coordinates
[tDesignTheta, tDesignPhi, ~] ...
    = cart2sph2(tDesign(:, 1), tDesign(:, 2), tDesign(:, 3));
[R2phi, R2rho] = KondorBackProj(tDesignTheta, tDesignPhi, scalingParam);
[R2x, R2y] = pol2cart2(R2phi, R2rho);

% Isolating t-design points falling within interval^2
tDesignInCube = R2x>=interval(1) & R2x<=interval(2) ...
    & R2y>=interval(1) & R2y<=interval(2);
R2x = R2x(tDesignInCube);
R2y = R2y(tDesignInCube);

% Create a grid in the cube interval^2
X = linspace(interval(1), interval(2), size(im, 1));
Y = linspace(interval(1), interval(2), size(im, 1));
[X, Y] = meshgrid(X, Y);

% Value for extrapulation outside interval^2
extval = 0;

% Compute the spherical harmonics matrix
if nargin<6
sh = sphericalHarmonics(tDesignTheta(tDesignInCube), ...
    tDesignPhi(tDesignInCube), ...
    bandlimit);
end

L = size(im, 3);
if L==1
    % Interpolate the image to values of the t-design that are back-projected
    % into the cube interval^2
    vals = interp2(X, Y, im, R2x, R2y, 'cubic', extval);

    % Estimate the spherical harmonics coefficients
    shc = (4*pi/size(tDesign, 1)) * (conj(sh) * vals);
else
    % Interpolate the image to values of the t-design that are back-projected
    % into the cube interval^2
    vals = zeros(length(R2x), L);
    parfor J=1:L
        vals(:, J) = interp2(X, Y, im(:, :, J), R2x, R2y, 'cubic', extval);
    end
    
    % Estimate the spherical harmonics coefficients
    shc = (4*pi/size(tDesign, 1)) * (conj(sh) * vals);
end
