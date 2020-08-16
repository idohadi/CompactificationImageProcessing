function rSHC = cSHC2rSHC(cSHC)
%%
% Call format
%   rSHC = cSHC2rSHC(cSHC)
% 
% Convert a complex column array containing spherical harmonics coefficents 
% to a real column array.
% 
%   rSHC(2*k-1) = real(cSHC(k))       and         rSHC(2*k) = imag(cSHC(k))
%   
% 
% Input arguments
%   cSHC            double      n x 1 complex array.
% 
% Output arguments
%   rSHC            double      n x 1 real array
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
rSHC = [real(cSHC), imag(cSHC)]';
rSHC = rSHC(:);
