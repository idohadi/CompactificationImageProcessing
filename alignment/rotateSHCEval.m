function rotatedSHC = rotateSHCEval(shc, N, rotation)
%%
% Call format
%   rotatedSHC = rotateSHCEval(shc, N, rotation)
%   values = fsht(shc, N, assocLegMat, phi, nthreads)
% 
% Calculate the spherical harmonics coefficients of the composition of a
% spherical function with a rotation.
% 
% This function rotated the spherical harmonics coefficients using
% the following procedure:
%   * Apply rotation to the meshgrid defined by theta and phi.
%   * Apply the forward spherical harmonics transform to compute the values
%     of the spherical function on the rotated grid.
%   * Calculate the spherical harmonics expansion of the rotated spherical 
%     function by applying the inverse spherical harmonics transform to the
%     computed values.
% 
% Input arguments
%   shc             double      (2*N)^2 x 1 array, 
%                               spherical harmonics coefficients.
%   N               double      positive integer, the bandlimit of shc is
%                               2*N-1.
%   rotation        double      4 x 1 array, a rotation in quaternion
%                               representation.

%   assocLegMat     double      cell array of matrices such that
%                               assocLegMat{m+2*N} is a (2*M-1) x (4*N-m)
%                               array satisfying
%                                   assocLegMat{m+2*N}(i, l-m+1) 
%                                       = >P_{l}^{|m|} (cos(theta(i)))
%                                   for 
%                                       0<=m<=2N-1 and m<=l<=2*N-1,
%                               where >P_{l}^{|m|} is a scaled, normalized
%                               associated Legendre function.
%   phi             double      K x 1 array, grid on [0,2pi].
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
% TODO
