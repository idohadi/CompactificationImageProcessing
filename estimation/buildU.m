function U = buildU(bandlimit, imageSize, tDesign, interval, scalingParam)
%%
% Call format
%   U = buildU(bandlimit, imageSize, tDesign, interval, scalingParam)
% 
% Compute the vector U which satisfies
%                           4*pi                          ( t-design points      )
%   U = SHC(image) = ---------------- * sphericalHarmonics( back-projected to    ) * P
%                    size(tDesign, 2)                     ( the cube interval^2  )
% where P is the interpolation operator computed by interp2linop.
% That is, U is the operator sending an image of size imageSize x imageSize
% to its spherical harmonics coefficients of order up to bandlimit.
% 
% Input arguments
%   bandlimit       double      positive integer, the bandlimit of the 
%                               sampled SHCs.
%   imageSize       double      positive integer, the image is of size 
%                                   imageSize x imageSize.
%   tDesign         double      N x 3 array, a spherical design (a 
%                               t-design) in Cartesian coordiantes.
%   interval        double      1 x 2 array, interval(1) is the lower bound
%                               of the interval and interval(2).
%   scalingParam    double      positive number, scaling parameter for 
%                               the projection.
%   
% Output arguments
%   U               double      (bandlimit+1)^2 x (imageSize^2) sparse 
%                               array, the matrix described above.
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

%% Compute U
% t-design in spherical coordinates
[tDesignTheta, tDesignPhi, ~] ...
    = cart2sph2(tDesign(:, 1), tDesign(:, 2), tDesign(:, 3));
[R2phi, R2rho] = KondorBackProj(tDesignTheta, tDesignPhi, scalingParam);
[R2x, R2y] = pol2cart2(R2phi, R2rho);

% Isolating t-design points falling within interval^2
tDesignInCube = R2x>=interval(1) & R2x<=interval(2) ...
    & R2y>=interval(1) & R2y<=interval(2);

% Compute the spherical harmonics matrix
sh = sphericalHarmonics(tDesignTheta(tDesignInCube), ...
    tDesignPhi(tDesignInCube), ...
    bandlimit);

% Compute U itself
U = (4*pi/size(tDesign, 1))*conj(sh)*interp2linop(imageSize, tDesign, interval, scalingParam);
U = sparse(U);
