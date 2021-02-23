function bispEst = mergeEstimatedBatches(batchBispEst, sampleSizes)
%%
% Call format
%   bispEst = mergeEstimatedBataches(batchBispEst, sampleSizes)
% 
% Given two estimators of the bispectrum, based on disjointed sample 
% sets,calculate the estimator based on the union of all samples.
% 
% Input arguments
%   batchBispEst        double      bispecturmLength x batchesNo array,
%                                   batchBispEst(:, j) is the bispectrum
%                                   estimator for the j-th batch.
%   sampleSizes         double      1 x batchesNo or batchesNo x 1 array, 
%                                   sampleSizes(j) is the sample size of
%                                   the j-th batch.
% 
% Output arguments
%   bispEst             double      k x 1 array, the estimator of the 
%                                   bispectrum based on the union of all 
%                                   sample sets.
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

%% Input validation
assert(size(batchBispEst, 2)==numel(sampleSizes), ...
    'Number of batches mismatch.');

%% Merge estimators
weights = sampleSizes(:);
weights = weights/sum(weights);
weights = weights';
bispEst = sum(batchBispEst.*weights, 2);
