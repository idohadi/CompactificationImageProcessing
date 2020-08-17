function [sphTheta, sphPhi] = KondorProj(R2phi, R2rho, scalingParam)
%%
% Call format
%   [sphTheta, sphPhi] = KondorProj(R2phi, R2rho, scalingParam)
% 
% Project points in R^2 in polar coordinates onto the 2-sphere parallel to 
% the z-axis. Output is in in spherical coordinates.
% 
% Convention
%   (x,y) is a point in R^2.
% 
%   Polar coordiantes:
%       Name    Interval
%       rho     [0, inifinity)
%       phi     [0, 2pi)
%   
%   (x,y), rho and phi satisfy
%       x = rho*cos(phi)
%       y = rho*sin(phi)
%   
%   (x,y,z) is a point in R^3.
%   Spherical coordiantes:
%       Variable    Interval
%       theta       [0,pi]
%       phi         [0,2pi)
%       rho         [0,infinity)
%   (x,y,z), theta, phi, rho satisfy
%       x = rho*cos(phi)*sin(theta)
%       y = rho*sin(phi)*sin(theta)
%       z = rho*cos(theta)
% 
% Input arguments
%   R2phi, R2rho        double      n x 1 arrays, points in R^2 in polar
%                                   coordiantes.
%   scalingParam        double      positive number, scaling parameter for 
%                                   the projection.
% 
% Output arguments
%   sphTheta, sphPhi    double      n x 1 arrays, points in S^2 in
%                                   spherical coordiantes.
% 
% Notes
%   (1) This function projects data from the plane onto the sphere,
%       parallel to the z-axis.
%   (2) For an example where this projection was used, see [1].
% 
% Reference
%   [1] Kakarala, R., & Mao, D. (2010). A theory of phase-sensitive 
%       rotation invariance with spherical harmonic and moment-based 
%       representations. 2010 IEEE Computer Society Conference on Computer 
%       Vision and Pattern Recognition, 105–112. 
%       https://doi.org/10.1109/CVPR.2010.5540222
% ***********************************************************
% Author    Ido Hadi
% Email     idohadi@mail.tau.ac.il
% Year      2020
% ***********************************************************
%% Project from R^2 to the sphere
sphTheta = scalingParam*R2rho;
sphPhi = R2phi;
