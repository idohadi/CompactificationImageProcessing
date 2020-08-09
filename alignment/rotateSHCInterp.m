function rotatedSHC = rotateSHCInterp(shc, N, rotation, assocLegMat, omega, theta, phi)
%%
% Call format
%   rotatedSHC = rotateSHCEval(shc, N, rotation)
%   values = fsht(shc, N, assocLegMat, phi, nthreads)
% 
% Calculate the spherical harmonics coefficients of the composition of a
% sphericla function with a rotation.
% 
% This function rotated the spherical harmonics coefficients using
% the following procedure:
%   * Apply rotation to the meshgrid defined by theta and phi.
%   * Evaluate the spherical function on the original grid. These are 
%     treated as its values on the rotated meshgrid.
%   * Interpolate these values to the meshgrid defined by theta and phi.
%   * Calculate the spherical harmonics expansion of the rotated spherical 
%     function by applying the inverse spherical harmonics transform to the
%     interpolated values.
% 
% Input arguments
%   shc             double      (2*N)^2 x 1 array, 
%                               spherical harmonics coefficients.
%   N               double      positive integer, the bandlimit of shc is
%                               2*N-1.
%   rotation        double      4 x 1 array, a rotation in quaternion
%                               representation.
%   assocLegMat     double      cell array of matrices as defined above.
%   omega           double      2N x 1 array, Guass-Legendre cubature
%                               weights.
%   theta           double      2N x 1 array, arccos of Guass-Legendre 
%                               cubature points.
%   phi             double      (4N-1) x 1 array, grid on [0, 2pi].
% 
% Output arguments
%   rotated SHC     double      spherical harmonics coefficients of the
%                               composition of the function expanded by shc
%                               with the rotation represented by rotation.
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
