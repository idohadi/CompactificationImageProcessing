function rotatedSHC = rotateSHC(shc, bandlimit, rotation, tDesign)
%%
% Call format
%   rotatedSHC = rotateSHC(shc, bandlimit, rotation, tDesign)
% 
% Calculate the spherical harmonics coefficients of the composition of a
% spherical function with a rotation.
% 
% This function rotated the spherical harmonics coefficients using
% the following procedure:
%   1. Apply rotation to the t-design.
%   2. Apply the forward spherical harmonics transform to compute the values
%      of the spherical function on the rotated t-design.
%      These values are the values of the rotated spherical function on the
%      original t-design.
%   3. Calculate the spherical harmonics expansion of the rotated spherical 
%      function by applying the inverse spherical harmonics transform to the
%      computed values.
% 
% Input arguments
%   shc         double      (bandlimit+1)^2 x 1 complex array, spherical 
%                           harmonics coefficients.
%   bandlimit   double      positive integer, the bandlimit of shc.
%   rotation    double      4 x 1 array, rotation in quaternion
%                           representation.
%   tDesign     double      N x 3 array, a spherical design (a t-design) in
%                           Cartesian coordiantes.
% 
% Output arguments
%   rotatedSHC  double      spherical harmonics coefficients of the
%                           composition of the function expanded by shc
%                           with the rotation represented by rotation.
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

%% Rotate SHC
% Step 1
if length(rotation)==4
    rotatedPoints = applyRotation(rotation, tDesign.');
elseif size(rotation, 1)==3 && size(rotation, 2)==3
    rotatedPoints = rotation*tDesign.';
else
    error('Rotation in unknown format.');
end
[theta, phi, ~] = cart2sph2(rotatedPoints(1, :), ...
                                            rotatedPoints(2, :), ...
                                            rotatedPoints(3, :));

% Step 2
sh = sphericalHarmonics(theta, phi, bandlimit);
rotatedFuncVals = sh.'*shc;

% Step 3
[theta, phi, ~] = cart2sph2(tDesign(:, 1), ...
                            tDesign(:, 2), ...
                            tDesign(:, 3));
sh = sphericalHarmonics(theta, phi, bandlimit);
rotatedSHC = (4*pi/size(tDesign, 1)) * (conj(sh) * rotatedFuncVals);
