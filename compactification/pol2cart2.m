function [x, y] = pol2cart2(phi, rho)
%%
% Call format
%   [x, y] = pol2cart2(phi, rho)
% 
% Convert the R^2 points from polar coordinates to Cartesian coordinates.
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
% Input arguments
%   phi, rho    double      n x 1 arrays, points in R^2 in polar 
%                           coordiantes.
% 
% Output arguments
%   x, y        double      n x 1 arrays, points in R^2 in Cartesian
%                           coordiantes.
% 
% Notes
%   This function depends on atan22 function from the repository
%   SmallRotationToolbox.
% 
% Reference
%   None
% ***********************************************************
% Author    Ido Hadi
% Email     idohadi@mail.tau.ac.il
% Year      2020
% ***********************************************************
%% Cartesian->polar
x = rho.*cos(phi);
y = rho.*sin(phi);
