function R = contractionMap(translation, rotation, scalingParam)
%%
% Call format
%   R = contractionMap(translation, rotation, scalingParam)
% 
% Applies the contraction map to an element of SE(2) to yield its
% corresponding element in SO(3).
% 
% Input arguments
%   translation     double      2 x 1 or 1 x 2 array, the translation.
%   rotation        double      scalar, the rotation angle in radians.
%   scalingParam    double      positive number, scaling parameter for the 
%                               projection.
% 
% Output arguments
%   R               double      3 x 3 array, the resulting element of SO(3)
%                               in matrix representation.
% 
% Notes
%   None
% 
% Reference
%   None
% ***********************************************************
% Author    Ido Hadi
% Email     idohadi@mail.tau.ac.il
% Year      2021
% ***********************************************************

%% Calculate contraction map
T = zeros(3, 3);
T(end, 1:2) = -translation;
T(1:2, end) = translation;
R = [cos(rotation), sin(rotation), 0; ...
    -sin(rotation), cos(rotation), 0; ...
    0,                  0,                  1];
R = expm(T/scalingParam)*R;
