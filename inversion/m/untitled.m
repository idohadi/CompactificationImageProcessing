%%
L = 3;
shc_bisp = calculateBispectrumOfRealValuedFunction(2*rand((L+1)^2, 1)-1, L);

problem.M = euclideanfactory((L+1)^2, 1);
problem.cost = @(x) sum(abs(calculateBispectrumOfRealValuedFunction(x, L)-shc_bisp))^2;
problem.grad = @(x) 2*calculateGradientOfBispectrumOfRealValuedFunction(x, L)*(calculateBispectrumOfRealValuedFunction(x, L)-shc_bisp);

checkgradient(problem);

%%
L = 3;
input = 2*rand((L+1)^2, 1)-1;

%%
for ind=23:23
    for k=1:(L+1)^2
        problem.M = euclideanfactory(1, 1);
problem.cost = @(x) testF(x, input, k, ind, L);
problem.grad = @(x) testG(x, input, k, ind, L);

checkgradient(problem);
disp(newline);
    end
end

%%
ind = 12;
k = 10;

problem.M = euclideanfactory(1, 1);
problem.cost = @(x) testF(x, input, k, ind, L);
problem.grad = @(x) testG(x, input, k, ind, L);

checkgradient(problem);
c = newline;

%%
for J=1:100
f = testF(rand(), input, k, ind, L)
g = testG(rand(), input, k, ind, L)
end


%%
L = 5; 
shc = 2*rand((L+1)^2, 1)-1;
shc_complexified = cfy(realSHC2ComplexSHC(shc, L));
shc_complexified = nshc(shc_complexified);
shc = complexSHC2RealSHC(rlfy(shc_complexified), L);

shc_bisp = calculateBispectrum(shc, L);


%%
[invertedSHC, squaredResidual, output] = invertRealValuedBispectrum(shc_bisp, L);

%%
shc_complexified = cfy(realSHC2ComplexSHC(shc, L));
invertedSHC_complexified = cfy(realSHC2ComplexSHC(invertedSHC, L));
alignedSHC2 = invertedSHC_complexified;

%%
[relativeDistance, alignedSHC2, optimalRotation, ...
    refinedMaxCorrelation, crudeMaxCorrelation, fminserachOutput] ...
    = alignSphericalHarmonics(shc_complexified, invertedSHC_complexified, L, td, 'SequenceSize', 2^8*72);

%%
pointNo = 10^5;
shcFunHandle = @(theta, phi) real(evalYt(theta, phi, L).' * nshc(shc_complexified));
recoveredFunHandle ...
    = @(theta, phi) real(evalYt(theta, phi, L).' * alignedSHC2);
rawRecoveredFunHandle ...
    = @(theta, phi) real(evalYt(theta, phi, L).' * invertedSHC_complexified);

figure;
tiledlayout(1, 3, ...
        'TileSpacing', 'normal', ...
        'Padding', 'Compact'); 

% Plot original estimated SHCs
nexttile;
plotSphericalFunction(shcFunHandle, pointNo);
title('Original SHC');
axis off;

% Plot recovered SHCs
nexttile;
plotSphericalFunction(recoveredFunHandle, pointNo);
title('Aligned recovered SHCs');
axis off;


nexttile;
plotSphericalFunction(rawRecoveredFunHandle , pointNo);
title('Raw recovered SHCs');
axis off;

% nexttile;
% plotSphericalFunction(originalRecoveredSHCFunHandle, pointNo);
% title('Misaligned recovered SHC');

%%
function f = testF(shc, input, k, ind, L)
input(k:k) = shc;
f = calculateBispectrumOfRealValuedFunction(input, L);
f = f(ind);
end


function grad = testG(shc, input, k, ind, L)
input(k:k) = shc;
grad = calculateGradientOfBispectrumOfRealValuedFunction(input, L);
grad = grad(k:k, ind);

% l1 = doesn't work;
%       m=11, 12 - don't work, the rest do
% l2 = works
% l  = works
end
