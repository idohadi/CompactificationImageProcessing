function [sampledImages, sampleVariances, basicViewOut] ...
    = generateImages(filename, imageSize, rotationAngles, translations, basicView)
%%
% Call format
%   generateImages(filename, imageSize, rotationAngles, translations)
%   generateImages(filename, imageSize, rotationAngles, translations, basicView)
%   sampledImages = generateImages(__)
%   [sampledImages, sampleVariances] = generateImages(__)
%   [sampledImages, sampleVariances, basicViewOut] = generateImages(__)
% 
% Generate a set of images, each is a Radon transform along the direction
% defined by basicView. The i-th image is rotated by rotationAngles(i) and
% translated by translation(:, i).
% 
% Input arguments
%   filename            double      char, path to .mrc file containing a
%                                   full 3d density map.
%   imageSize           double      positive integer, the resulting images
%                                   would be imageSize x imageSize arrays.
%   rotationAngles      double      n-elements array, a list of angles in
%                                   radians.
%   translations        double      2 x n array, a list of shifts in
%                                   pixels.
%   basicView           double      3 x 3 array, a rotation applied to the
%                                   map such that the images are taken 
%                                   w.r.t. to the rotated map.
% 
% Output arguments
%   sampledImages       double      imageSize x imageSize x n array, the
%                                   resulting image array.
%   sampleVariances     double      n x 1 array, the variance of the
%                                   correspodning images.
%   basicViewOut        double      3 x 3 array, a rotation matrix.
% 
% Notes
%   This function depends on the packages ASPIRE and SmallRotationToolbox.
% 
% Reference
%   None
% ***********************************************************
% Author    Ido Hadi
% Email     idohadi@mail.tau.ac.il
% Year      2020
% ***********************************************************

%% Handle input
narginchk(4, 5);

% Handle basicView optional argument
if nargin==4
    basicView = quaternion2rotMat(randRotation(1));
end
basicViewOut = basicView;

assert(ischar(filename) | isstring(filename), 'File name must be a char or string.');
assert(isscalar(imageSize) & isnumeric(imageSize) ...
    & imageSize>0 & round(imageSize)==imageSize, ...
    'Image size must be a positive integer scalar.');
assert(isnumeric(rotationAngles), 'Rotation angles must be numeric.');
assert(isnumeric(translations) & size(translations, 1)==2, ...
    'Translation must be numeric 2 x n array.');
assert(isnumeric(basicView) & all(size(basicView)==[3, 3]), ...
    'Basic view must be a 3 x 3 rotation matrix.');

assert(size(translations, 2)==numel(rotationAngles), ...
    'There must be equal number of translations and rotations.');

rotationAngles = rotationAngles(:);
n = length(rotationAngles); % Sample size

%% Generate images
% Load the density map
map = ReadMRC(filename);

% Convert rotations to rotation matrices
rotations = zeros([3, 3, n]);
rotations(3, 3, :) = 1;
rotations(1, 1, :) = cos(rotationAngles);
rotations(2, 2, :) = cos(rotationAngles);
rotations(2, 1, :) = -sin(rotationAngles);
rotations(1, 2, :) = sin(rotationAngles);
% Incorporate the basic view
rotations = reshape(rotations, [3, 3*n]);
rotations = basicView*rotations;
rotations = reshape(rotations, [3, 3, n]);
rotations = permute(rotations, [2, 1, 3]);

% Generate the sampled images
sampledImages = cryo_project(map, rotations, imageSize);
sampledImages = cryo_addshifts(sampledImages, translations');

% Calculate the variances
sampleVariances = squeeze(var(sampledImages, 0, [1, 2]));
