%%
addpath 'C:\Code\SphericalInvariantAnalysis\bin\mex'
addpath 'C:\Code\SphericalInvariantAnalysis\inversion\m'
addpath 'C:\Code\SphericalInvariantAnalysis\spherical_spectra\m'
addpath 'C:\Code\SphericalInvariantAnalysis\spherical_harmonics\m'

%% Test: complex-valued function inversion
% Test parameters
Ls = [4, 6, 8, 10, 12];
td = loadtd('sf025.00339');
trialsNo = 20;

% Test results struct
results = struct();

% The test
for ell=1:length(Ls)
    for trial=1:trialsNo
        shc = cfy(randomNormalizedSHC(Ls(ell), 1));
        shcBisp = calculateBispectrum_M(shc, Ls(ell));

        t = tic();
        [invertedSHC, rootedResidual, output] = invertComplexValuedBispectrum(shcBisp, Ls(ell));
        inversionRuntime = toc(t);
        
        t = tic();
        [relativeDistance, alignedSHC2, optimalRotation, refinedMaxCorrelation, crudeMaxCorrelation, fminserachOutput] = alignSphericalHarmonics(shc, invertedSHC, Ls(ell), td);
        alignmentRuntime = toc(t);

        % Save results
        results(trial, ell).bandlimit = Ls(ell);
        results(trial, ell).shc = shc;
        results(trial, ell).shcBisp = shcBisp;
        results(trial, ell).invertedSHC = invertedSHC;
        results(trial, ell).rootedResidual = rootedResidual;
        results(trial, ell).output = output;
        results(trial, ell).relativeDistance = relativeDistance;
        results(trial, ell).alignedSHC2 = alignedSHC2;
        results(trial, ell).optimalRotation = optimalRotation;
        results(trial, ell).refinedMaxCorrelation = refinedMaxCorrelation;
        results(trial, ell).crudeMaxCorrelation = crudeMaxCorrelation;
        results(trial, ell).fminserachOutput = fminserachOutput;
        results(trial, ell).inversionRuntime = inversionRuntime;
        results(trial, ell).alignmentRuntime = alignmentRuntime;

        disp(['iter #', num2str(trial), ' of ', num2str(trialsNo), ...
            '. Bandlimit = ', num2str(Ls(ell)), ...
            '. rel dist = ', num2str(relativeDistance), ...
            '. Inv runtime = ', num2str(inversionRuntime), ...
            '. Align runtime = ', num2str(alignmentRuntime), '.']);
    end
end
complexValuedSHCInversionTestResultsNewMEXFunc = results;
save('complexValuedSHCInversionTestResultsNewMEXFunc.mat', 'complexValuedSHCInversionTestResultsNewMEXFunc');

%% Test: real-valued function inversion 
results = load('complexValuedSHCInversionTestResultsNewMEXFunc.mat');
results = results.complexValuedSHCInversionTestResultsNewMEXFunc;
Ls = [4, 6, 8, 10, 12];
td = loadtd('sf025.00339');

resultsFollowup = cell(length(Ls), 1);
for J=1:size(results, 2)
    res = results([results(:, J).relativeDistance]>=0.1, J);
    resultsFollowup{J} = struct();
    
    for K=1:length(res)
        t = tic();
        rootedResidual = Inf;
        while rootedResidual>10^-2
            [invertedSHC, rootedResidual, output] = invertComplexValuedBispectrum(res(K).shcBisp, res(K).bandlimit);
        end
        inversionRuntime = toc(t);
        
        t = tic();
        [relativeDistance, alignedSHC2, optimalRotation, refinedMaxCorrelation, crudeMaxCorrelation, fminserachOutput] = alignSphericalHarmonics(res(K).shc, invertedSHC, res(K).bandlimit, td);
        alignmentRuntime = toc(t);

        % Save results
        resultsFollowup{J}(K).bandlimit = res(K).bandlimit;
        resultsFollowup{J}(K).shc = res(K).shc;
        resultsFollowup{J}(K).shcBisp = res(K).shcBisp;
        resultsFollowup{J}(K).invertedSHC = res(K).invertedSHC;
        resultsFollowup{J}(K).rootedResidual = rootedResidual;
        resultsFollowup{J}(K).output = output;
        resultsFollowup{J}(K).relativeDistance = relativeDistance;
        resultsFollowup{J}(K).alignedSHC2 = alignedSHC2;
        resultsFollowup{J}(K).optimalRotation = optimalRotation;
        resultsFollowup{J}(K).refinedMaxCorrelation = refinedMaxCorrelation;
        resultsFollowup{J}(K).crudeMaxCorrelation = crudeMaxCorrelation;
        resultsFollowup{J}(K).fminserachOutput = fminserachOutput;
        resultsFollowup{J}(K).inversionRuntime = inversionRuntime;
        resultsFollowup{J}(K).alignmentRuntime = alignmentRuntime;

        disp(['iter #', num2str(K), ' of ', num2str(length(res)), ...
            '. Bandlimit = ', num2str(res(K).bandlimit), ...
            '. rel dist = ', num2str(relativeDistance), ...
            '. Inv runtime = ', num2str(inversionRuntime), ...
            '. Align runtime = ', num2str(alignmentRuntime), '.']);
    end
end

complexValuedSHCInversionTestResultsNewMEXFuncFollowup = results;
save('complexValuedSHCInversionTestResultsNewMEXFuncFollowup.mat', 'complexValuedSHCInversionTestResultsNewMEXFuncFollowup');


%% Test: real-valued function inversion
% Test parameters
Ls = [4, 6, 8, 10, 12];
td = loadtd('sf025.00339');
trialsNo = 20;

% Test results struct
results = struct();

% The test
for ell=1:length(Ls)
    for trial=1:trialsNo
        shc = r2c(randomNormalizedSHC(Ls(ell), 0), Ls(ell));
        shcBisp = calculateBispectrum_M(shc, Ls(ell));

        t = tic();
        [invertedSHC, rootedResidual, output] = invertComplexValuedBispectrum(shcBisp, Ls(ell));
        inversionRuntime = toc(t);

        t = tic();
        [relativeDistance, alignedSHC2, optimalRotation, refinedMaxCorrelation, crudeMaxCorrelation, fminserachOutput] = alignSphericalHarmonics(shc, invertedSHC, Ls(ell), td);
        alignmentRuntime = toc(t);

        % Save results
        results(trial, ell).bandlimit = Ls(ell);
        results(trial, ell).shc = shc;
        results(trial, ell).shcBisp = shcBisp;
        results(trial, ell).invertedSHC = invertedSHC;
        results(trial, ell).rootedResidual = rootedResidual;
        results(trial, ell).output = output;
        results(trial, ell).relativeDistance = relativeDistance;
        results(trial, ell).alignedSHC2 = alignedSHC2;
        results(trial, ell).optimalRotation = optimalRotation;
        results(trial, ell).refinedMaxCorrelation = refinedMaxCorrelation;
        results(trial, ell).crudeMaxCorrelation = crudeMaxCorrelation;
        results(trial, ell).fminserachOutput = fminserachOutput;
        results(trial, ell).inversionRuntime = inversionRuntime;
        results(trial, ell).alignmentRuntime = alignmentRuntime;
        
        disp(['iter #', num2str(trial), ' of ', num2str(trialsNo), ...
            '. Bandlimit = ', num2str(Ls(ell)), ...
            '. rel dist = ', num2str(relativeDistance), ...
            '. Inv runtime = ', num2str(inversionRuntime), ...
            '. Align runtime = ', num2str(alignmentRuntime), '.']);
    end
end
realValuedSHCInversionTestResultsNewMEXFunc = results;
save('realValuedSHCInversionTestResultsNewMEXFunc.mat', 'realValuedSHCInversionTestResultsNewMEXFunc');
