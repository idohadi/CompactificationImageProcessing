function cSHC = rSHC2cSHC(rSHC)
%%
% Call format
%   cSHC = rSHC2cSHC(rSHC)
% 
% Convert a real column array containing spherical harmonics coefficents to
% a compelx column array.
% 
%   rSHC(2*k-1) = real(cSHC(k))       and         rSHC(2*k) = imag(cSHC(k))
%   
% 
% Input arguments
%   rSHC            double      n x 1 real array
% 
% Output arguments
%   cSHC            double      n x 1 complex array.
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
%% Convert
cSHC = reshape(rSHC, [2, round(length(rSHC)/2)])';
cSHC = complex(cSHC(:, 1) + 1i*cSHC(:, 2));
