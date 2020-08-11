function CGvec = ClebschGordanCoeffs(l1, l2, l, m)
%%
% Call format
%   shc = randSHC(bandlimit)
%   shc = randSHC(bandlimit, N)
% 
% Calculate the Clebsch-Gordan coefficeints of order (l1, l2, l, m), where
%   l1>=0, l2>=, abs(l1-l2)<=l<=l1+l2, abs(m)<=l.
% 
% Input arguments
%   l1, l2, l, m        double      integers, representing the order of the
%                                   Clebsch-Gordan coefficients.
% 
% Output arguments
%   CGvec               double      N x 1 array, such that 
%                                   CGvec(m1 - m1min + 1) is the 
%                                   Clebsch-Gordan coefficient 
%                                       <l1 m1 l2 (m-m1) | l m >
%                                   where 
%                                       m1min = max{-l1, m-l2}, 
%                                       m1max = min{l1, m+l2}, and
%                                       N = m1max - m1min + 1.
% 
% Notes
%   (1) This function depends on ClebschGordanCoeffs_mex.
%   (2) The coefficients are calcualted using  the method of [1]. In [1]
%       there is also a concise description of theoretical background.
% 
% Reference
%   [1] Straub, W. O. (n.d.). Efficient Computation of Clebsch-Gordan 
%       Coefficients. Retrieved October 28, 2019, 
%       from http://vixra.org/abs/1403.0263
% ***********************************************************
% Author    Ido Hadi
% Email     idohadi@mail.tau.ac.il
% Year      2020
% ***********************************************************
%% Input validation
assert(l1>=0 && l2>=0, 'l1 and l2 must be non-negative.');
assert(abs(l1-l2)<=l && l<=l1+l2, ...
    'l1, l2 and l must satisfy abs(l1-l2)<=l<=l1+l2.');
assert(abs(m)<=l, 'l and m must satisfy abs(m)<=l.');

%% Calculation
CGvec = ClebschGordanCoeffs_mex(l1, l2, l, m);
