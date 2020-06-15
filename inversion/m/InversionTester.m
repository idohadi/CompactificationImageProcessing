%% Parameters
L = 8;
td = loadtd('sf071.02593');

%%
L = 6; 
td = loadtd('sf071.02593');
REAL=0;
shc = randomNormalizedSHC(L, REAL); ans
[b, g] = calculateBispectrum(shc, L, REAL); 

shcRot = c2r(rotateSphericalHarmonicsByEstimation(r2c(shc, L), L, td, generateUniformlyRandomRotations(1)), L);
t = tic();
norm(calculateBispectrum(shc, L) - calculateBispectrum(shcRot, L))/norm(calculateBispectrum(shc, L))
toc(t)

A = calculateBispectrum(shc, L); 
A = A(2:2:end);

%%
K = bisp_vec(r2c(shc, L), L, CGs);
K = imag(K);

%%
[invertedSHC, rootedResidual, output] = invertRealValuedBispectrum(b, L);

%%
shcComplexified = r2c(shc, L);
invertedSHCComplexified = r2c(invertedSHC, L);

%%
[relativeDistance, alignedSHC2, optimalRotation, refinedMaxCorrelation, crudeMaxCorrelation, fminserachOutput] = alignSphericalHarmonics(shcComplexified, invertedSHCComplexified, L, td, 'SequenceSize', 2^9*72);

%%
pointNo = 10^4;
imFunHandle = @(theta, phi) real(evalYt(theta, phi, L).' * shcComplexified);
recoveredFunHandle ...
    = @(theta, phi) real(evalYt(theta, phi, L).' * invertedSHCComplexified);

figure;
tiledlayout(1, 2, ...
        'TileSpacing', 'normal', ...
        'Padding', 'Compact');

% Plot original SHCs
nexttile;
caxis([-3, 3]);
plotSphericalFunction(imFunHandle, pointNo);
title('Estimated SHCs of image');
axis off;

% Plot recovered SHCs
nexttile;
caxis([-3, 3]);
plotSphericalFunction(recoveredFunHandle, pointNo);
title('Recovered SHCs');
axis off;

%% Real-valued SHC alignment test
trialsNo = 100;

results = struct();
for trial=1:trialsNo
    shc1 = 2*rand((L+1)^2, 1)-1;
    shc1Complexified = r2c(shc1, L);
    shc1Complexified = nshc(shc1Complexified);
    shc1 = c2r(shc1Complexified, L);

    trueRotation = generateUniformlyRandomRotations(1);
    shc2Complexified ...
        = rotateSphericalHarmonicsByEstimation(shc1Complexified, L, td, trueRotation);
    shc2 = c2r(shc2Complexified, L);
    
    t = tic();
    [relativeDistance, alignedSHC2Complexified, optimalRotation, ...
    refinedMaxCorrelation, crudeMaxCorrelation, fminserachOutput] ...
        = alignSphericalHarmonics(shc1Complexified, shc2Complexified, ...
            L, td, 'SequenceSize', 2^9*72);
    alignedSHC2 = c2r(alignedSHC2Complexified, L);
    runtime = toc(t);
    
    % Save results
    results(trial).shc1 = shc1;
    results(trial).shc2 = shc2;
    results(trial).shc1Complexified = shc1Complexified;
    results(trial).shc2Complexified = shc2Complexified;
    results(trial).trueRotation = trueRotation;

    results(trial).preAlignmentRelativeDistance ...
        = norm(shc1Complexified-shc2Complexified)/norm(shc1Complexified);
    results(trial).postAlignmentRelativeDistance = relativeDistance;
    results(trial).alignedSHC2 = alignedSHC2;
    results(trial).alignedSHC2Complexified = alignedSHC2Complexified;
    results(trial).trueRotation = optimalRotation;
    results(trial).refinedMaxCorrelation = refinedMaxCorrelation;
    results(trial).crudeMaxCorrelation = crudeMaxCorrelation;
    results(trial).fminserachOutput = fminserachOutput;
    results(trial).runtime = runtime;
    
    disp(['Trial #', num2str(trial), '. ', ...
        'Pre-alignment rel. dist. = ', ...
            num2str(norm(shc1Complexified-shc2Complexified)/norm(shc1Complexified)), ', ', ...
        'Post-alignment rel. dist = ', num2str(relativeDistance), ...
        '. Runtime (sec) = ', num2str(runtime), '.']);
end

realValuedSHCAlignmentTestResults = results;
save('realValuedSHCAlignmentTestResults.mat', 'realValuedSHCAlignmentTestResults');


%% Parameters
L = 10;
td = loadtd('sf071.02593');

opts = optimoptions(@lsqnonlin, ...
    'SpecifyObjectiveGradient', true, ...
    'OptimalityTolerance', 10^-12, ...
    'FunctionTolerance', 10^-15, ...
    'StepTolerance', 10^-10, ...
    'Display', 'off'); 

%% Real-valued SHC bispectrum inversion test
experimentsNo = 5;
trialsNo = 75;

results = struct();

for experiment=1:experimentsNo
    shc = 2*rand((L+1)^2, 1)-1;
    shcComplexified = r2c(shc, L);
    shcComplexified = nshc(shcComplexified);
    shc = c2r(shcComplexified, L);

    shcBispectrum = calculateBispectrumOfRealValuedFunction(shc, L);
    
    for trial=1:trialsNo
        x0 = rand((L+1)^2, 1);
        x0Complexified = r2c(x0, L);
        x0Complexified = nshc(x0Complexified);
        x0 = c2r(x0Complexified, L);
        func = @(shc) inversionObjectiveFunc(shc, shcBispectrum, L);
        
        t = tic();
        [invertedSHC, res, ~, ~, output] = lsqnonlin(func, x0, [], [], opts);
        bispectrumInversionRuntime = toc(t);
        rootedResidual = sqrt(res);
        
        invertedSHCComplexified = r2c(invertedSHC, L);
        
        t = tic();
        [relativeDistance, alignedInvertedSHCComplexified, optimalRotation, ...
            refinedMaxCorrelation, crudeMaxCorrelation, fminserachOutput] ...
            = alignSphericalHarmonics(shcComplexified, invertedSHCComplexified, ...
                    L, td, 'SequenceSize', 2^9*72);
        alignmentRuntime = toc(t);
        
        alignedInvertedSHC = c2r(alignedInvertedSHCComplexified, L);
        
        % Save results
        results(trial, experiment).shc = shc;
        results(trial, experiment).shcComplexified = shcComplexified;
        results(trial, experiment).shcBispectrum = shcBispectrum;
        results(trial, experiment).x0 = x0;
        results(trial, experiment).rootedResidual = rootedResidual;
        results(trial, experiment).invertedSHC = invertedSHC;
        results(trial, experiment).invertedSHCComplexified ...
            = invertedSHCComplexified;
        results(trial, experiment).output = output;
        results(trial, experiment).bispectrumInversionRuntime ...
            = bispectrumInversionRuntime;
        results(trial, experiment).alignedInvertedSHCComplexified ...
            = alignedInvertedSHCComplexified;
        results(trial, experiment).alignedInvertedSHC = alignedInvertedSHC;
        results(trial, experiment).alignmentRuntime = alignmentRuntime;
        
        disp(['Exp. #', num2str(experiment), ...
        ', Trial #', num2str(trial), '. ', ...
        'Pre-alignment rel. dist. = ', ...
            num2str(norm(shcComplexified-invertedSHCComplexified)/norm(shcComplexified)), ', ', ...
        'Post-alignment rel. dist = ', num2str(relativeDistance), ...
        '. Buisp. Inv. Runtime (sec) = ', num2str(bispectrumInversionRuntime), '.', ...
        'Align. Runtime (sec) = ', num2str(alignmentRuntime), ]);
    end
end

realValuedBispectrumInversionTestResults = results;
save('realValuedBispectrumInversionTestResults.mat', 'realValuedBispectrumInversionTestResults');


%%
% Alignment test
markerSize = 3;
trialNo = length(realValuedSHCAlignmentTestResults);

fig1 = figure();

t = tiledlayout(2, 1, 'Padding', 'compact', 'TileSpacing', 'compact');

nexttile;
hold;
preAlign = [realValuedSHCAlignmentTestResults.preAlignmentRelativeDistance];
postAlign = [realValuedSHCAlignmentTestResults.postAlignmentRelativeDistance];
plot(1:trialNo, preAlign, ...
    'Marker', 'o', ...
    'MarkerFaceColor', [0 0.4470 0.7410], ...
    'MarkerSize', markerSize);
plot(1:trialNo, postAlign, ...
    'Marker', 'o', ...
    'MarkerFaceColor', [0.8500 0.3250 0.0980], ...
    'MarkerSize', markerSize);
hold off;

title('Aligning Randomly Rotated SHC, Real-Valued Function');
xlabel('Trial');
ylabel('Relative error');
legend({'Pre-alignment', 'Post-alignment'});
set(gca, 'yscale', 'log');

nexttile;
runtime = [realValuedSHCAlignmentTestResults.runtime];
plot(1:trialNo, runtime, ...
    'Marker', 'o', ...
    'MarkerFaceColor', [0 0.4470 0.7410], ...
    'MarkerSize', markerSize);

% title('Aligning Randomly Rotated SHC, Real-Valued Function');
xlabel('Trial');
ylabel('Runtime (sec)');


% Bispectrum inversion, relative error
[trialNo, experimentNo] = size(realValuedBispectrumInversionTestResults);
colors = [0 0.4470 0.7410; ...
            0.8500 0.3250 0.0980; ...
            0.9290 0.6940 0.1250; ...
            0.4940 0.1840 0.5560; ...
            0.4660 0.6740 0.1880];
fig2 = figure();

t = tiledlayout(5, 1, 'Padding', 'compact', 'TileSpacing', 'compact');
title(t, {'Bispectrum Inversion, Real-Valued Function,', 'Relative Error'});

relativeError = zeros(trialNo, experimentNo);
for J=1:experimentNo
    for N=1:trialNo
        alignedInvertedSHC ...
            = realValuedBispectrumInversionTestResults(N, J).alignedInvertedSHC;
        shc = realValuedBispectrumInversionTestResults(N, J).shc;
        
        relativeError(N, J) = norm(alignedInvertedSHC - shc)/norm(shc);
    end
end

for J=1:experimentNo
    nexttile;
    plot(1:trialNo, relativeError(:, J), ...
        'Color', colors(J, :), ...
        'Marker', 'o', ...
        'MarkerFaceColor', colors(J, :), ...
        'MarkerSize', markerSize);
    
    xlim([1, trialNo]);
    xticks(0:5:trialNo);
    ylim([0.6, 1.1]);
    yticks(0.6:0.1:1.1);
    xlabel('Trial');
%     ylabel('Relative error');
    title(['Experiment #', num2str(J)]);
end

% Bispectrum inversion, non-squared residual
fig3 = figure();
t = tiledlayout(5, 1, 'Padding', 'compact', 'TileSpacing', 'compact');
title(t, {'Bispectrum Inversion, Real-Valued Function,' 'Non-squared Residual'});

residual = zeros(trialNo, experimentNo);
for J=1:experimentNo
    for N=1:trialNo
        residual(N, J) ...
            = realValuedBispectrumInversionTestResults(N, J).rootedResidual;
    end
end

for J=1:experimentNo
    nexttile;
    plot(1:trialNo, residual(:, J), ...
        'Color', colors(J, :), ...
        'Marker', 'o', ...
        'MarkerFaceColor', colors(J, :), ...
        'MarkerSize', markerSize);
    set(gca, 'yscale', 'log');
    
    xlim([1, trialNo]);
    xticks(0:5:trialNo);
    
    ylim([10^-10, 1]);
    yticks(10.^(-10:1));
    
    xlabel('Trial');
%     ylabel('Non-squared residual');
    title(['Experiment #', num2str(J)]);
end

% Bispectrum inversion, first order optimality condition
fig4 = figure();
t = tiledlayout(5, 1, 'Padding', 'compact', 'TileSpacing', 'compact');
title(t, {'Bispectrum Inversion, Real-valued Function,', 'First Order Optimality Condition'});

firstOrderOptimality = zeros(trialNo, experimentNo);
for J=1:experimentNo
    for N=1:trialNo
        output = realValuedBispectrumInversionTestResults(N, J).output;
        firstOrderOptimality(N, J) = output.firstorderopt;
    end
end

for J=1:experimentNo
    nexttile;
    plot(1:trialNo, firstOrderOptimality(:, J), ...
        'Color', colors(J, :), ...
        'Marker', 'o', ...
        'MarkerFaceColor', colors(J, :), ...
        'MarkerSize', markerSize);
    set(gca, 'yscale', 'log');
    
    xlim([1, trialNo]);
    xticks(0:5:trialNo);
    
    ylim([10^-11, 10^-4]);
    yticks(10.^(-10:-4));
    
    xlabel('Trial');
%     ylabel('First order opt. cond.');
    title(['Experiment #', num2str(J)]);
end

% Bispectrum inversion, Runtime
fig5 = figure();
t = tiledlayout(5, 1, 'Padding', 'compact', 'TileSpacing', 'compact');
title(t, {'Bispectrum Inversion,', 'Runtime (sec)'});

bispectrumInversionRuntime = zeros(trialNo, experimentNo);
alignmentRuntime = zeros(trialNo, experimentNo);
for J=1:experimentNo
    for N=1:trialNo
        bispectrumInversionRuntime(N, J) ... 
            = realValuedBispectrumInversionTestResults(N, J).bispectrumInversionRuntime;
        alignmentRuntime(N, J) ... 
            = realValuedBispectrumInversionTestResults(N, J).alignmentRuntime;
    end
end

for J=1:experimentNo
    nexttile;
    
    xlabel('Trial');
    xlim([1, trialNo]);
    xticks(0:5:trialNo);
    title(['Experiment #', num2str(J)]);
    
    yyaxis left;
    plot(1:trialNo, bispectrumInversionRuntime(:, J), ...
        'Color', colors(1, :), ...
        'Marker', 'o', ...
        'MarkerFaceColor', colors(1, :), ...
        'MarkerSize', markerSize);
    ylabel('Bispectrm inversion');
%     ylim([10^-11, 10^-4]);
%     yticks(10.^(-10:-4));
    
    yyaxis right;
    plot(1:trialNo, alignmentRuntime(:, J), ...
        'Color', colors(2, :), ...
        'Marker', 'o', ...
        'MarkerFaceColor', colors(2, :), ...
        'MarkerSize', markerSize);
    ylabel('Alignment');
%     ylim([10^-11, 10^-4]);
%     yticks(10.^(-10:-4));
    
end


%% Utility functions
function complexCoeffs = r2c(realCoeffs, bandlimit)
complexCoeffs = zeros(2*(bandlimit+1)^2, 1);

for order=0:bandlimit
    for degree=0:order
        if degree==0
            complexCoeffs(2*(order^2+order)+1) ...
                = realCoeffs(order^2+1);
        else
            % Handling +degree
            % Real part of coefficient
            complexCoeffs(2*(order^2+order+degree)+1) ...
                = realCoeffs(order^2+1+2*(degree-1)+1);
            % Imaginary part of coefficeint
            complexCoeffs(2*(order^2+order+degree)+2) ...
                = realCoeffs(order^2+1+2*(degree-1)+2);
            
            % Handling -degree
            % Real part of coefficient
            complexCoeffs(2*(order^2+order-degree)+1) ...
                = (-1)^degree * realCoeffs(order^2+1+2*(degree-1)+1);
            % Imaginary part of coefficeint
            complexCoeffs(2*(order^2+order-degree)+2) ...
                = (-1)^(degree+1) * realCoeffs(order^2+1+2*(degree-1)+2);
        end
    end
end

complexCoeffs = reshape(complexCoeffs, [2, (bandlimit+1)^2]).';
complexCoeffs = complexCoeffs(:, 1) + 1i*complexCoeffs(:, 2);

end

function realCoeffs = c2r(complexCoeffs, bandlimit)
complexCoeffs = reshape([real(complexCoeffs), imag(complexCoeffs)].', [2*length(complexCoeffs), 1]);

realCoeffs = zeros((bandlimit+1)^2, 1);
for order=0:bandlimit
    for degree=0:order
        if degree==0
            realCoeffs(order^2+1) ...
                = complexCoeffs(2*(order^2+order)+1);
        else
            % Handling +degree
            % Real part of coefficient
            realCoeffs(order^2+1+2*(degree-1)+1) ...
                = complexCoeffs(2*(order^2+order+degree)+1);
            
            % Imaginary part of coefficeint
             realCoeffs(order^2+1+2*(degree-1)+2) ...
                 = complexCoeffs(2*(order^2+order+degree)+2);
            
        end
    end
end
end

function [F, grad] = inversionObjectiveFunc(shc, bispectrum, L)
F = calculateBispectrumOfRealValuedFunction(shc, L) - bispectrum;
grad = calculateGradientOfBispectrumOfRealValuedFunction(shc, L)';
end