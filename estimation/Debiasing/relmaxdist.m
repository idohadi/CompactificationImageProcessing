function reldist = relmaxdist(ests, truth)
% Measures the distance between vector ests to truth in the max-norm
% relative to truth.
% 
% Input arguments:
%   ests    -   (L+1)^2 vector; estimator of SHC
%   truth   -   (L+1)^2 vector; the real SHC
% 
% Output arguments:
%   reldist -   the max-norm of ests-truth dided by the max-norm of truth
% 

reldist = max(abs(ests - truth), [], 1)/max(abs(truth));
