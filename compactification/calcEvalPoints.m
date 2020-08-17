function evalPoints = calcEvalPoints(theta, phi, scalingParam)
%%
% Call format
%   evalPoints = calcEvalPoints(theta, phi)
% 
% Given a quadrature on the sphere, creates a meshgrid with theta on the
% y-axis and phi on the x-axis, back-project it onto R^2, and then reshapes
% it to be linear in a column-major order.
% 
% Input arguments
%   theta           double      2N x 1 array, arccos of Guass-Legendre 
%                               cubature points.
%   phi             double      (4N-1) x 1 array, grid on [0, 2pi].
%   scalingParam    double      positive number, scaling parameter for the 
%                               projection.
% 
% Output arguments
%   evalPoints      double      
% 
% Notes
%   None
% 
% Reference
%   None
% ***********************************************************
% Author    Ido Hadi
% Email     idohadi@mail.tau.ac.il
% Year      2020
% ***********************************************************
%% Compute evaluation points
[phi, theta] = meshgrid(phi, theta);
[R2phi, R2rho] = KondorBackProj(theta, phi, scalingParam);
[x, y] = pol2cart2(R2phi, R2rho);
evalPoints = [x(:), y(:)];
