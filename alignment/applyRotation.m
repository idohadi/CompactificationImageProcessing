function xRot = applyRotation(rotation, x)
%%
% Apply the rotation to vectors in R^3.
% 
% Input arguments
%   rotation    double      4 x 1 or 1 x 4 real array, a unit quaternion
%                           representing the rotation.
%   x           double      3 x N or N x 3 real array.
% 
% Output arguments
%   xRot        double      3 x N or N x 3 real array.
% 
% Notes
%   (1) This is a wrapper of the mex function applyRotation_mex.
% 
% Reference
%   None
% ***********************************************************
% Author    Ido Hadi
% Email     idohadi@mail.tau.ac.il
% Year      2020
% ***********************************************************
%% Input validation
assert(ismatrix(x), 'x must be a matrix.');
assert(~isempty(x), 'x can''t be empty matrix.');
if size(x, 1)~=3
    assert(size(x, 2)==3, 'x must have a dimension of size 3');
    x = x';
end

%% Calling the mex function
applyRotation_mex(rotation, x);
