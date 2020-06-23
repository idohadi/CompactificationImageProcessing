%% Parameters
L = 8;
td = loadtd('sf020.00222');
CGs = build_CGs_vec(L);

%% Test: complex-valued function inversion
trialsNo = 20;
results = struct();

for trial=1:trialsNo
    shc = cfy(randomNormalizedSHC(L, 1));
    
    shcBisp = bisp_vec(shc, L, CGs);
    
    [f, res, output] = invbisp(shcBisp, L, CGs);
    [relativeDistance, alignedSHC2, optimalRotation, refinedMaxCorrelation, crudeMaxCorrelation, fminserachOutput] = alignSphericalHarmonics(shc, f, L, td);
    
    % Save results
    results(trial).bandlimit = L;
    results(trial).shc = shc;
    results(trial).shcBisp = shcBisp;
    results(trial).invertedSHC = f;
    results(trial).rootedResidual = sqrt(res);
    results(trial).output = output;
    results(trial).relativeDistance = relativeDistance;
    results(trial).alignedSHC2 = alignedSHC2;
    results(trial).optimalRotation = optimalRotation;
    results(trial).refinedMaxCorrelation = refinedMaxCorrelation;
    results(trial).crudeMaxCorrelation = crudeMaxCorrelation;
    results(trial).fminserachOutput = fminserachOutput;
    
    disp(['iter #', num2str(trial), ' of ', num2str(trialsNo), ...
        '. rel dist = ', num2str(relativeDistance), '.']);
end

complexValuedSHCInversionTestResultsNewTD = results;
save('complexValuedSHCInversionTestResultsNewTD.mat', 'complexValuedSHCInversionTestResultsNewTD');

%% Test: complex-valued function inversion, effect of initial point
L = 8;
td = loadtd('sf020.00222');
CGs = build_CGs_vec(L);

expNo = 5;
trialsNo = 500;
results = struct();

for expInd=1:expNo
    shc = cfy(randomNormalizedSHC(L, 1));
    
    shcBisp = bisp_vec(shc, L, CGs);
    b = shcBisp;
    
    for trial=1:trialsNo
        opts = optimoptions(@lsqnonlin, ...
        'SpecifyObjectiveGradient', true, ...
        'Display', 'iter', ...
        'CheckGradients', true); 

        func = @(x) objFunc(cfy(x), L, b, CGs);

        % Generate random initial point for optimization.
        % It is normalized to have powspec of all ones
        x0 = randomNormalizedSHC(L, 1);

        % Run optimization algorithm
        t = tic();
        [f, res, ~, ~, output] = lsqnonlin(func, x0, [], [], opts);
        inversionRuntime = toc(t);
        f = cfy(f);
        
        % Align SHCs
        t = tic();
        [relativeDistance, alignedSHC2, optimalRotation, refinedMaxCorrelation, crudeMaxCorrelation, fminserachOutput] = alignSphericalHarmonics(shc, f, L, td);
        alignRuntime = toc(t);
        
        % Save things
        results(trial, expInd).opts = opts;
        
        results(trial, expInd).shc = shc;
        results(trial, expInd).shcBisp = shcBisp;
        results(trial, expInd).L = L;
        results(trial, expInd).x0 = x0;
        results(trial, expInd).rootedResidual = sqrt(res);
        results(trial, expInd).output = output;
        results(trial, expInd).opts = opts;
        
        results(trial, expInd).relativeDistance = relativeDistance;
        results(trial, expInd).alignedSHC2 = alignedSHC2;
        results(trial, expInd).optimalRotation = optimalRotation;
        results(trial, expInd).refinedMaxCorrelation = refinedMaxCorrelation;
        results(trial, expInd).crudeMaxCorrelation = crudeMaxCorrelation;
        results(trial, expInd).fminserachOutput = fminserachOutput;
        
        results(trial, expInd).alignRuntime = alignRuntime;
        results(trial, expInd).inversionRuntime = inversionRuntime;
        
        disp(['iter #', num2str(trial), ' of ', num2str(trialsNo), ...
            '. rel dist = ', num2str(relativeDistance), '. ', ...
            'Inv Runtime = ', num2str(inversionRuntime), '. ', ...
            'Align Runtime = ', num2str(alignRuntime), '.']);
    end
end

complexValuedSHCInverstionInitialPointTestNewTD = results;
save('complexValuedSHCInverstionInitialPointTestNewTD.mat', 'complexValuedSHCInverstionInitialPointTestNewTD');