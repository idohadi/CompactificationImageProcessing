%%
addpath 'C:\Code\SphericalInvariantAnalysis\bin\mex'
addpath 'C:\Code\SphericalInvariantAnalysis\inversion\m'
addpath 'C:\Code\SphericalInvariantAnalysis\spherical_spectra\m'
addpath 'C:\Code\SphericalInvariantAnalysis\spherical_harmonics\m'

%% Test: complex-valued function inversion
% Test parameters
Ls = [5, 6, 8, 10, 12];
td = loadtd('sf025.00339');
trialsNo = 20;

% Test results struct
results = struct();

% The test
for L=Ls
    for trial=1:trialsNo
        shc = cfy(randomNormalizedSHC(L, 1));
        shcBisp = calculateBispectrum_M(shc, L);

        t = tic();
        [invertedSHC, rootedResidual, output] = invertComplexValuedBispectrum(shcBisp, L);
        inversionRuntime = toc(t);

        t = tic();
        [relativeDistance, alignedSHC2, optimalRotation, refinedMaxCorrelation, crudeMaxCorrelation, fminserachOutput] = alignSphericalHarmonics(shc, invertedSHC, L, td);
        alignmentRuntime = toc(t);

        % Save results
        results(trial).bandlimit = L;
        results(trial).shc = shc;
        results(trial).shcBisp = shcBisp;
        results(trial).invertedSHC = invertedSHC;
        results(trial).rootedResidual = rootedResidual;
        results(trial).output = output;
        results(trial).relativeDistance = relativeDistance;
        results(trial).alignedSHC2 = alignedSHC2;
        results(trial).optimalRotation = optimalRotation;
        results(trial).refinedMaxCorrelation = refinedMaxCorrelation;
        results(trial).crudeMaxCorrelation = crudeMaxCorrelation;
        results(trial).fminserachOutput = fminserachOutput;
        results(trial).inversionRuntime = inversionRuntime;
        results(trial).alignmentRuntime = alignmentRuntime;

        disp(['iter #', num2str(trial), ' of ', num2str(trialsNo), ...
            '. Bandlimit = ', num2str(L), ...
            '. rel dist = ', num2str(relativeDistance), ...
            '. Inv runtime = ', num2str(inversionRuntime), ...
            '. Align runtime = ', num2str(alignmentRuntime), '.']);
    end
end
complexValuedSHCInversionTestResultsNewMEXFunc = results;
save('complexValuedSHCInversionTestResultsNewMEXFunc.mat', 'complexValuedSHCInversionTestResultsNewMEXFunc');

%% Test: real-valued function inversion
% Test parameters
Ls = [5, 6, 8, 10, 12];
td = loadtd('sf025.00339');
trialsNo = 20;

% Test results struct
results = struct();

% The test
for L=Ls
    for trial=1:trialsNo
        shc = cfy(r2c(randomNormalizedSHC(L, 0)));
        shcBisp = calculateBispectrum_M(shc, L);

        t = tic();
        [invertedSHC, rootedResidual, output] = invertComplexValuedBispectrum(shcBisp, L);
        inversionRuntime = toc(t);

        t = tic();
        [relativeDistance, alignedSHC2, optimalRotation, refinedMaxCorrelation, crudeMaxCorrelation, fminserachOutput] = alignSphericalHarmonics(shc, invertedSHC, L, td);
        alignmentRuntime = toc(t);

        % Save results
        results(trial).bandlimit = L;
        results(trial).shc = shc;
        results(trial).shcBisp = shcBisp;
        results(trial).invertedSHC = invertedSHC;
        results(trial).rootedResidual = rootedResidual;
        results(trial).output = output;
        results(trial).relativeDistance = relativeDistance;
        results(trial).alignedSHC2 = alignedSHC2;
        results(trial).optimalRotation = optimalRotation;
        results(trial).refinedMaxCorrelation = refinedMaxCorrelation;
        results(trial).crudeMaxCorrelation = crudeMaxCorrelation;
        results(trial).fminserachOutput = fminserachOutput;
        results(trial).inversionRuntime = inversionRuntime;
        results(trial).alignmentRuntime = alignmentRuntime;

        disp(['iter #', num2str(trial), ' of ', num2str(trialsNo), ...
            '. rel dist = ', num2str(relativeDistance), ...
            '. Inv runtime = ', num2str(inversionRuntime), ...
            '. Align runtime = ', num2str(alignmentRuntime), '.']);
    end
end
realValuedSHCInversionTestResultsNewMEXFunc = results;
save('realValuedSHCInversionTestResultsNewMEXFunc.mat', 'realValuedSHCInversionTestResultsNewMEXFunc');
