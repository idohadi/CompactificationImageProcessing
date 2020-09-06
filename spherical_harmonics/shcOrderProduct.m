function factoredSHC = shcOrderProduct(shc, bandlimit, factors)
%%
% Call format
%   factoredSHC = shcOrderProduct(shc, bandlimit, factors)
% 
% Return spherical harmonics coefficeints such that the coefficents of
% order l are multiplied by factrors(l+1). That is,
%   out(l*(l+1)+m+1) = factors(l+1)*shc(l*(l+1)+m+1);
% 
% 
% Input arguments
%   shc             double      (bandlimit+1)^2 x N complex array, 
%                               spherical harmonics coefficients.
%   bandlimit       double      positive integer, the bandlimit of shc.
%   factors         double      (bandlimit+1) x 1 complex array.
% 
% Output arguments
%   factoredSHC     double      the output described above.
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

%% Factor product
factoredSHC = zeros(size(shc));
for l=0:bandlimit
    for m=-l:l
        factoredSHC(l*(l+1)+m+1) = factors(l+1)*shc(l*(l+1)+m+1);
    end
end
