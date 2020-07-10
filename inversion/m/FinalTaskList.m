%% Final Task List
addpath 'C:\Code\SphericalInvariantAnalysis\bin\mex'
addpath 'C:\Code\SphericalInvariantAnalysis\inversion\m'
addpath 'C:\Code\SphericalInvariantAnalysis\spherical_spectra\m'
addpath 'C:\Code\SphericalInvariantAnalysis\spherical_harmonics\m'

%% Task 2
% In the real-valued case, what is the effect of the distance between the 
% initial point and the ground truth?
%% 
L = 6;
td = loadtd('sf025.00339');
trialsNo = 100;
expNo = 10;
pow_lb = -5;
pow_ub = 0;

opts = optimoptions(@lsqnonlin, ...
    'SpecifyObjectiveGradient', true, ...
    'OptimalityTolerance', 10^-12, ...
    'FunctionTolerance', 10^-15, ...
    'StepTolerance', 10^-10, ...
    'Display', 'off'); 
tol = 10^-8;

% Test results struct
results = struct();

% The test
for ex=1:expNo
    shc = r2c(randomNormalizedSHC(L, 0), L);
    shcBisp = calculateBispectrum_M(shc, L);
    
    for trial=1:trialsNo
        t = tic();
        func = @(x) inversionObjectiveFunc(cfy(x), shcBisp, L);
        % Generating initial point
        direction = r2c(randomNormalizedSHC(L, 0), L);
        direction = direction/norm(direction, 2);
        digit = randi(9);
        pow = randi(6) - 6;
        initDistance = digit*10^pow;
        x0 = rlfy(initDistance*direction + shc);
        % Inverting the bipsectrm of shc
        [invertedSHC, rootedResidual, ~, ~, output] = lsqnonlin(func, x0, [], [], opts);
        invertedSHC = cfy(invertedSHC);
        rootedResidual = sqrt(rootedResidual);
        inversionRuntime = toc(t);
        
        
        % Aligning the inverted SHC and the original one
        t = tic();
        [relativeDistance, alignedSHC2, optimalRotation, refinedMaxCorrelation, crudeMaxCorrelation, fminserachOutput] = alignSphericalHarmonics(shc, invertedSHC, L, td);
        alignmentRuntime = toc(t);

        % Save results
        results(trial, ex).bandlimit = L;
        results(trial, ex).shc = shc;
        results(trial, ex).shcBisp = shcBisp;
        results(trial, ex).invertedSHC = invertedSHC;
        results(trial, ex).rootedResidual = rootedResidual;
        results(trial, ex).x0 = x0;        
        results(trial, ex).output = output;
        results(trial, ex).relativeDistance = relativeDistance;
        results(trial, ex).alignedSHC2 = alignedSHC2;
        results(trial, ex).optimalRotation = optimalRotation;
        results(trial, ex).refinedMaxCorrelation = refinedMaxCorrelation;
        results(trial, ex).crudeMaxCorrelation = crudeMaxCorrelation;
        results(trial, ex).fminserachOutput = fminserachOutput;
        results(trial, ex).inversionRuntime = inversionRuntime;
        results(trial, ex).alignmentRuntime = alignmentRuntime;
        results(trial, ex).initDistance = initDistance;
        results(trial, ex).direction = direction;

        disp(['exp #', num2str(ex), ' of ', num2str(expNo), ...
            '. iter #', num2str(trial), ' of ', num2str(trialsNo), ...
            '. Bandlimit = ', num2str(L), ...
            '. Inv runtime = ', num2str(inversionRuntime), ...
            '. Align runtime = ', num2str(alignmentRuntime), '.']);
        disp(['. init dist = ', num2str(initDistance), ...
            '. final res = ', num2str(rootedResidual), ...
            '. rel dist = ', num2str(relativeDistance)]);
    end
end

Task2L6 = results;
save('Task2L6.mat', 'Task2L6');
%%
% ===================
% Ploting the results
% ===================

results = load('Task2L6.mat');
results = results.Task2L6;

color =     [0, 0.4470, 0.7410; ...
            [0.8500 0.3250 0.0980]; ...
            [0.9290 0.6940 0.1250]; ...
            [0.4940 0.1840 0.5560]; ...
            [0.4660 0.6740 0.1880]];
sz = 3;

fig = figure;
tiledlayout(3, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile;
hold;
for J=1:expNo
    x = [results(:, J).initDistance];
    x = x(:);
    [x, inds] = sortrows(x);
    y = [results(:, J).relativeDistance];
    y = y(:);
    y = y(inds);
    plot(x, y, ...
        'Marker', 'o', ...
        'MarkerSize', sz);
end
hold off;
set(gca, 'xscale', 'log');
set(gca, 'yscale', 'log');

title('Bispectrum inversion of real-valued spherical function');
% xlabel('Distance between initial guess and ground truth');
ylabel('Relative error');
yticks(10.^(-13:1));

% ============

nexttile;
hold;
for J=1:expNo
    S = [results(:, J).output];
    
    x = [results(:, J).initDistance];
    x = x(:);
    [x, inds] = sortrows(x);
    y = [S.firstorderopt];
    y = y(:);
    y = y(inds);
    plot(x, y, ...
        'Marker', 'o', ...
        'MarkerSize', sz);
end
hold off;

set(gca, 'yscale', 'log');
set(gca, 'xscale', 'log');

% xlabel('Distance between initial guess and ground truth');
yticks(10.^(-13:1));
ylabel('Residual (rooted)');

% ============

nexttile;

hold;
for J=1:expNo
    S = [results(:, J).output];
    
    x = [results(:, J).initDistance];
    x = x(:);
    [x, inds] = sortrows(x);
    y = [S.firstorderopt];
    y = y(:);
    y = y(inds);
    plot(x, y, ...
        'Marker', 'o', ...
        'MarkerSize', sz);
end
hold off;
set(gca, 'xscale', 'log');
set(gca, 'yscale', 'log');

% xlabel('Distance between initial guess and ground truth');
ylabel('First order optimality condition');
yticks(10.^(-13:1));

legcell = cell(1, size(results, 2));
for J=1:size(results, 2)
    legcell{J} = ['Experiment #', num2str(J)];
end
legend(legcell, ...
    'Location', 'southoutside', ...
    'Orientation', 'horizontal', ...
    'NumColumns', 5);

save2pdf('Task2L6fig.pdf', fig);

%% Task 4
% Add a small random imaginary perturbation to the real-valued spherical 
% function. Can one invert the bispectrum of that?
%% 
L = 6;
td = loadtd('sf025.00339');
trialsNo = 100;

opts = optimoptions(@lsqnonlin, ...
    'SpecifyObjectiveGradient', true, ...
    'OptimalityTolerance', 10^-12, ...
    'FunctionTolerance', 10^-15, ...
    'StepTolerance', 10^-10, ...
    'Display', 'off'); 
tol = 10^-8;

% Test results struct
results = struct();

% The test
for trial=1:trialsNo
    shc = r2c(randomNormalizedSHC(L, 0), L);
    shcBisp = calculateBispectrum_M(shc, L);
    
    t = tic();
    % Pertrub shc
    direction = zeros((L+1)^2, 1);
    direction(1) = (2*randi(2)-3)*1i;
    digit = randi(9);
    pow = randi(6) - 6;
    initDistance = digit*10^pow;
    shcPertrubed = initDistance*direction + shc;
    shcPertrubedBisp = calculateBispectrum_M(shcPertrubed, L);
    
    % Inverting the bipsectrm of shc
    x0 = randomNormalizedSHC(L, 1);
    func = @(x) inversionObjectiveFunc(cfy(x), shcPertrubedBisp, L);
    [invertedSHC, rootedResidual, ~, ~, output] = lsqnonlin(func, x0, [], [], opts);
    invertedSHC = cfy(invertedSHC);
    rootedResidual = sqrt(rootedResidual);
    inversionRuntime = toc(t);

    % Aligning the inverted SHC and the original one
    t = tic();
    [relativeDistance, alignedSHC2, optimalRotation, refinedMaxCorrelation, crudeMaxCorrelation, fminserachOutput] = alignSphericalHarmonics(shcPertrubed, invertedSHC, L, td);
    alignmentRuntime = toc(t);

    % Save results
    results(trial).bandlimit = L;
    results(trial).shc = shc;
    results(trial).shcPertrubed = shcPertrubed;
    results(trial).shcBisp = shcBisp;
    results(trial).shcPertrubedBisp = shcPertrubedBisp;
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
    results(trial).initDistance = initDistance;
    results(trial).direction = direction;

    disp(['iter #', num2str(trial), ' of ', num2str(trialsNo), ...
        '. Bandlimit = ', num2str(L), ...
        '. Inv runtime = ', num2str(inversionRuntime), ...
        '. Align runtime = ', num2str(alignmentRuntime), '.']);
    disp(['first ord. opt. = ', num2str(output.firstorderopt), ...
        '. rooted residual = ', num2str(rootedResidual), ...
        '. rel dist = ', num2str(relativeDistance), ...
        '. init dist = ', num2str(initDistance)]);
    disp('');
end

Task4L6 = results;
save('Task4L6.mat', 'Task4L6');

%%
% ===================
% Ploting the results
% ===================

results = load('Task4L6.mat');
results = results.Task4L6;

color =     [0, 0.4470, 0.7410; ...
            [0.8500 0.3250 0.0980]; ...
            [0.9290 0.6940 0.1250]; ...
            [0.4940 0.1840 0.5560]; ...
            [0.4660 0.6740 0.1880]];
sz = 3;

fig = figure;
tiledlayout(3, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile;
hold;
x = [results.initDistance];
x = x(:);
[x, inds] = sortrows(x);
y = [results.relativeDistance];
y = y(:);
y = y(inds);
plot(x, y, ...
    'Marker', 'o', ...
    'MarkerSize', sz);
hold off;
set(gca, 'xscale', 'log');

title('Bispectrum inversion of real-valued spherical function');
% xlabel('Pertubation size');
ylabel('Relative error');
yticks(0:0.1:1);

% ============

nexttile;
hold;
x = [results.initDistance];
x = x(:);
[x, inds] = sortrows(x);
y = [results.rootedResidual];
y = y(:);
y = y(inds);
plot(x, y, ...
    'Marker', 'o', ...
    'MarkerSize', sz);
hold off;
set(gca, 'xscale', 'log');
set(gca, 'yscale', 'log');

% xlabel('Pertubation size');
ylabel('Residual (rooted)');
yticks(10.^(-13:1));

% ============

nexttile;

hold;
S = [results.output];

x = [results.initDistance];
x = x(:);
[x, inds] = sortrows(x);
y = [S.firstorderopt];
y = y(:);
y = y(inds);
plot(x, y, ...
    'Marker', 'o', ...
    'MarkerSize', sz);
hold off;
set(gca, 'xscale', 'log');
set(gca, 'yscale', 'log');

xlabel('Pertubation size');
ylabel('First order optimality condition');
yticks(10.^(-13:1));

save2pdf('Task4L6fig.pdf', fig);


%% Task 5
% Add a small random imaginary perturbation to the real-valued spherical 
% function. Can one invert the bispectrum of that?
%% 
L = 6;
td = loadtd('sf025.00339');
trialsNo = 100;

opts = optimoptions(@lsqnonlin, ...
    'SpecifyObjectiveGradient', true, ...
    'OptimalityTolerance', 10^-12, ...
    'FunctionTolerance', 10^-15, ...
    'StepTolerance', 10^-10, ...
    'Display', 'off'); 
tol = 10^-8;
maxTries = 10;

% Test results struct
results = struct();

% The test
for trial=1:trialsNo
    shc = r2c(randomNormalizedSHC(L, 0), L);
    shcBisp = calculateBispectrum_M(shc, L);
    
    t = tic();
    % Pertrub shc
    direction = cfy(randomNormalizedSHC(L, 1));
    direction = direction/norm(direction);
    digit = randi(9);
    pow = randi(6) - 6;
    initDistance = digit*10^pow;
    shcPertrubed = initDistance*direction + shc;
    shcPertrubedBisp = calculateBispectrum_M(shcPertrubed, L);
    
    % Inverting the bipsectrm of shc
    r = Inf;
    tries = 0;
    bestResidual = Inf;
    bestSHC = shc;
    while tries<10
        x0 = randomNormalizedSHC(L, 1);
        func = @(x) inversionObjectiveFunc(cfy(x), shcPertrubedBisp, L);
        [invertedSHC, rootedResidual, ~, ~, output] = lsqnonlin(func, x0, [], [], opts);
        
        r = output.firstorderopt;
        tries = tries + 1;
        
        if bestResidual>sqrt(rootedResidual)
            bestRediaul = sqrt(rootedResidual);
            bestSHC = invertedSHC;
            bestOutput = output;
        end
        
        if sqrt(rootedResidual)<10^-8 && r<=tol
            break;
        end
    end
    invertedSHC = cfy(bestSHC);
    rootedResidual = sqrt(rootedResidual);
    inversionRuntime = toc(t);

    % Aligning the inverted SHC and the original one
    t = tic();
    [relativeDistance, alignedSHC2, optimalRotation, refinedMaxCorrelation, crudeMaxCorrelation, fminserachOutput] = alignSphericalHarmonics(shcPertrubed, invertedSHC, L, td);
    alignmentRuntime = toc(t);

    % Save results
    results(trial).bandlimit = L;
    results(trial).shc = shc;
    results(trial).shcPertrubed = shcPertrubed;
    results(trial).shcBisp = shcBisp;
    results(trial).shcPertrubedBisp = shcPertrubedBisp;
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
    results(trial).initDistance = initDistance;
    results(trial).direction = direction;
    results(trial).tries = tries;
    results(trial).maxTries = maxTries;
    results(trial).bestRediaul = bestRediaul;
    results(trial).bestSHC = bestSHC;
    results(trial).bestOutput = bestOutput;
    
    
    disp(['iter #', num2str(trial), ' of ', num2str(trialsNo), ...
        '. Bandlimit = ', num2str(L), ...
        '. Inv runtime = ', num2str(inversionRuntime), ...
        '. Align runtime = ', num2str(alignmentRuntime), '.']);
    disp(['first ord. opt. = ', num2str(output.firstorderopt), ...
        '. rooted residual = ', num2str(rootedResidual), ...
        '. rel dist = ', num2str(relativeDistance), ...
        '. init dist = ', num2str(initDistance), ...
        '. tries = ', num2str(tries)]);
    disp('');
end

Task5L6 = results;
save('Task5L6.mat', 'Task5L6');

%%
% ===================
% Ploting the results
% ===================

results = load('Task5L6.mat');
results = results.Task5L6;

color =     [0, 0.4470, 0.7410; ...
            [0.8500 0.3250 0.0980]; ...
            [0.9290 0.6940 0.1250]; ...
            [0.4940 0.1840 0.5560]; ...
            [0.4660 0.6740 0.1880]];
sz = 3;

fig = figure;
tiledlayout(3, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile;
hold;
x = [results.initDistance];
x = x(:)/sqrt(7);
[x, inds] = sortrows(x);
y = [results.relativeDistance];
y = y(:);
y = y(inds);
plot(x, y, ...
    'Marker', 'o', ...
    'MarkerSize', sz);
hold off;
set(gca, 'xscale', 'log');
set(gca, 'yscale', 'log');

title('Bispectrum inversion of real-valued spherical function');
% xlabel('Relative pertubation size');
ylabel('Relative error');
yticks(10.^(-13:1));

% ============

nexttile;
hold;
x = [results.initDistance];
x = x(:)/sqrt(7);
[x, inds] = sortrows(x);
y = [results.rootedResidual];
y = y(:);
y = y(inds);
plot(x, y, ...
    'Marker', 'o', ...
    'MarkerSize', sz);
hold off;
set(gca, 'xscale', 'log');
set(gca, 'yscale', 'log');

% xlabel('Relative pertubation size');
ylabel('Residual (rooted)');
yticks(10.^(-13:1));

% ============

nexttile;

hold;
S = [results.output];

x = [results.initDistance];
x = x(:)/sqrt(7);
[x, inds] = sortrows(x);
y = [S.firstorderopt];
y = y(:);
y = y(inds);
plot(x, y, ...
    'Marker', 'o', ...
    'MarkerSize', sz);
hold off;
set(gca, 'xscale', 'log');
set(gca, 'yscale', 'log');

xlabel('Relative pertubation size');
ylabel('First order optimality condition');
yticks(10.^(-13:1));

save2pdf('Task5L6fig.pdf', fig);


%% Task 3
% What happens when I use a MATLAB  function that takes into account the 
% symmetries in the real-valued case?
%% 
L = 6;
td = loadtd('sf025.00339');
trialsNo = 100;

opts = optimoptions(@lsqnonlin, ...
    'SpecifyObjectiveGradient', true, ...
    'OptimalityTolerance', 10^-12, ...
    'FunctionTolerance', 10^-15, ...
    'StepTolerance', 10^-10, ...
    'Display', 'off', ...
    'CheckGradients', true);
tol = 10^-8;
maxTries = 10;
[bispInds, symComp] = utilityFunc(L);

% Test results struct
results = struct();

% The test
for trial=1:trialsNo
    shc = r2c(randomNormalizedSHC(L, 0), L);
    shcBisp = calculateBispectrum_M(shc, L);
    
    % Inverting the bipsectrm of shc
    t = tic();
    r = Inf;
    tries = 0;
    bestResidual = Inf;
    bestSHC = shc;
    while tries<10
        x0 = randomNormalizedSHC(L, 0);
        func = @(x) inversionObjectiveFuncReal(r2c(x, L), shcBisp, L, bispInds, symComp);
        [invertedSHC, rootedResidual, ~, ~, output] = lsqnonlin(func, x0, [], [], opts);
        
        r = rootedResidual;
        tries = tries + 1;
        
        if bestResidual>sqrt(rootedResidual)
            bestRediaul = sqrt(rootedResidual);
            bestSHC = r2c(invertedSHC, L);
            bestOutput = output;
            initResidual = norm(func(x0), 2);
        end
        
        if sqrt(rootedResidual)<10^-8 && r<=tol
            break;
        end
    end
    invertedSHC = bestSHC;
    rootedResidual = bestRediaul;
    inversionRuntime = toc(t);
    
    % Aligning the inverted SHC and the original one
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
    results(trial).initResidual = initResidual;
    results(trial).optimalRotation = optimalRotation;
    results(trial).refinedMaxCorrelation = refinedMaxCorrelation;
    results(trial).crudeMaxCorrelation = crudeMaxCorrelation;
    results(trial).fminserachOutput = fminserachOutput;
    results(trial).inversionRuntime = inversionRuntime;
    results(trial).alignmentRuntime = alignmentRuntime;
    results(trial).tries = tries;
    results(trial).maxTries = maxTries;
    results(trial).bestRediaul = bestRediaul;
    results(trial).bestSHC = bestSHC;
    results(trial).bestOutput = bestOutput;
    
    
    disp(['iter #', num2str(trial), ' of ', num2str(trialsNo), ...
            '. Bandlimit = ', num2str(L), ...
            '. Inv runtime = ', num2str(inversionRuntime), ...
            '. Align runtime = ', num2str(alignmentRuntime), '.']);
    disp(['. init dist = ', num2str(initDistance), ...
        '. init res = ', num2str(initResidual), ...
        '. final res = ', num2str(rootedResidual), ...
        '. rel dist = ', num2str(relativeDistance)]);
end

Task3L6 = results;
save('Task3L6.mat', 'Task3L6');

%%
% ===================
% Ploting the results
% ===================

results = load('Task3L6.mat');
results = results.Task3L6;

color =     [0, 0.4470, 0.7410; ...
            [0.8500 0.3250 0.0980]; ...
            [0.9290 0.6940 0.1250]; ...
            [0.4940 0.1840 0.5560]; ...
            [0.4660 0.6740 0.1880]];
sz = 3;

fig = figure;
tiledlayout(3, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile;
hold;
y = [results.relativeDistance];
y = y(:);
plot(y, ...
    'Marker', 'o', ...
    'MarkerSize', sz);
hold off;
set(gca, 'yscale', 'log');

title('Bispectrum inversion of real-valued spherical function');
% xlabel('Relative pertubation size');
ylabel('Relative error');
yticks(0:0.1:1);

% ============

nexttile;
hold;
J = 1;
y = [results.rootedResidual];
y = y(:);
y2 = [results.initResidual];
y2 = y2(:);

yyaxis left;
plot(y, ...
    'Marker', 'o', ...
    'MarkerSize', sz);
set(gca, 'yscale', 'log');
ylabel('Final residual (rooted)');
yticks(10.^(-13:1));

yyaxis right;
plot(y2, ...
    'Marker', 'o', ...
    'MarkerSize', sz);
hold off;
yticks(0:1:10);
% xlabel('Relative pertubation size');
ylabel('Initial residual (rooted)');

% ============

nexttile;

hold;
S = [results.output];

y = [S.firstorderopt];
y = y(:);
plot(y, ...
    'Marker', 'o', ...
    'MarkerSize', sz);
hold off;
set(gca, 'yscale', 'log');

xlabel('Trial');
ylabel('First order optimality condition');
yticks(10.^(-13:1));

% save2pdf('Task3L6fig.pdf', fig);

%% Task 2, with the code from task 3
% In the real-valued case, what is the effect of the distance between the 
% initial point and the ground truth?
%% 
L = 6;
td = loadtd('sf025.00339');
trialsNo = 100;
expNo = 1;
pow_lb = -5;
pow_ub = 0;

opts = optimoptions(@lsqnonlin, ...
    'SpecifyObjectiveGradient', true, ...
    'OptimalityTolerance', 10^-12, ...
    'FunctionTolerance', 10^-15, ...
    'StepTolerance', 10^-10, ...
    'Display', 'off'); 
tol = 10^-8;
[bispInds, symComp] = utilityFunc(L);

% Test results struct
results = struct();

% The test
for ex=1:expNo
    shc = randomNormalizedSHC(L, 0);
    shcBisp = calculateBispectrum_M(r2c(shc, L), L);
    
    for trial=1:trialsNo
        t = tic();
        problem.cost = @(x) inversionObjectiveFunc(r2c(x, L), shcBisp, L, bispInds, symComp);
        % Generating initial point
        direction = randomNormalizedSHC(L, 0);
        direction = direction/norm(direction, 2);
        digit = randi(9);
        pow = randi(6) - 6;
        initDistance = digit*10^pow;
        x0 = initDistance*direction + shc;
        % Inverting the bipsectrm of shc
        [invertedSHC, rootedResidual, ~, ~, output] = lsqnonlin(func, x0, [], [], opts);
        invertedSHC = r2c(invertedSHC, L);
        rootedResidual = sqrt(rootedResidual);
        inversionRuntime = toc(t);
        
        % Aligning the inverted SHC and the original one
        t = tic();
        [relativeDistance, alignedSHC2, optimalRotation, refinedMaxCorrelation, crudeMaxCorrelation, fminserachOutput] = alignSphericalHarmonics(r2c(shc, L), invertedSHC, L, td);
        alignmentRuntime = toc(t);

        % Save results
        results(trial, ex).bandlimit = L;
        results(trial, ex).shc = r2c(shc, L);
        results(trial, ex).shcBisp = shcBisp;
        results(trial, ex).invertedSHC = invertedSHC;
        results(trial, ex).rootedResidual = rootedResidual;
        results(trial, ex).output = output;
        results(trial, ex).relativeDistance = relativeDistance;
        results(trial, ex).alignedSHC2 = alignedSHC2;
        results(trial, ex).optimalRotation = optimalRotation;
        results(trial, ex).refinedMaxCorrelation = refinedMaxCorrelation;
        results(trial, ex).crudeMaxCorrelation = crudeMaxCorrelation;
        results(trial, ex).fminserachOutput = fminserachOutput;
        results(trial, ex).inversionRuntime = inversionRuntime;
        results(trial, ex).alignmentRuntime = alignmentRuntime;
        results(trial, ex).initDistance = initDistance;
        results(trial, ex).direction = direction;

        disp(['exp #', num2str(ex), ' of ', num2str(expNo), ...
            '. iter #', num2str(trial), ' of ', num2str(trialsNo), ...
            '. Bandlimit = ', num2str(L), ...
            '. rel dist = ', num2str(relativeDistance), ...
            '. init dist = ', num2str(initDistance), ...
            '. Inv runtime = ', num2str(inversionRuntime), ...
            '. Align runtime = ', num2str(alignmentRuntime), '.']);
    end
end

Task2w3L6 = results;
save('Task2w3L6.mat', 'Task2w3L6');

%%
% ===================
% Ploting the results
% ===================

results = load('Task2w3L6.mat');
results = results.Task2w3L6;

color =     [0, 0.4470, 0.7410; ...
            [0.8500 0.3250 0.0980]; ...
            [0.9290 0.6940 0.1250]; ...
            [0.4940 0.1840 0.5560]; ...
            [0.4660 0.6740 0.1880]];
sz = 3;

fig = figure;
tiledlayout(3, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile;
hold;
for J=1:expNo
    x = [results(:, J).initDistance];
    x = x(:);
    [x, inds] = sortrows(x);
    y = [results(:, J).relativeDistance];
    y = y(:);
    y = y(inds);
    plot(x, y, ...
        'Marker', 'o', ...
        'MarkerSize', sz);
end
hold off;
set(gca, 'xscale', 'log');
set(gca, 'yscale', 'log');

title('Bispectrum inversion of real-valued spherical function');
% xlabel('Distance between initial guess and ground truth');
ylabel('Relative error');
yticks(10.^(-13:1));

% ============

nexttile;
hold;
for J=1:expNo
    x = [results(:, J).initDistance];
    x = x(:);
    [x, inds] = sortrows(x);
    y = [results(:, J).rootedResidual];
    y = y(:);
    y = y(inds);
    
    
    initialResidual = [];
    for K=1:size(x)
        initialResidual = norm(results(K, J).shcBisp, 2);
    end
    plot(x, y, ...
        'Marker', 'o', ...
        'MarkerSize', sz);
end
hold off;
set(gca, 'xscale', 'log');
set(gca, 'yscale', 'log');

% xlabel('Distance between initial guess and ground truth');
ylabel('Residual (rooted)');
yticks(10.^(-13:1));

% ============

nexttile;

hold;
for J=1:expNo
    S = [results(:, J).output];
    
    x = [results(:, J).initDistance];
    x = x(:);
    [x, inds] = sortrows(x);
    y = [S.firstorderopt];
    y = y(:);
    y = y(inds);
    plot(x, y, ...
        'Marker', 'o', ...
        'MarkerSize', sz);
end
hold off;
set(gca, 'xscale', 'log');
set(gca, 'yscale', 'log');

xlabel('Distance between initial guess and ground truth');
ylabel('First order optimality condition');
yticks(10.^(-13:1));

save2pdf('Task2w3L6fig.pdf', fig);


%% Task 7
% Try using Manopt as a solver of the optimization problem on the 
% appropriate product manifold
%% 
L = 6;
td = loadtd('sf025.00339');
trialsNo = 100;

opts = optimoptions(@lsqnonlin, ...
    'SpecifyObjectiveGradient', true, ...
    'OptimalityTolerance', 10^-12, ...
    'FunctionTolerance', 10^-15, ...
    'StepTolerance', 10^-10, ...
    'Display', 'off');
tol = 10^-8;
maxTries = 10;

% Manopt manifold
mans = struct();
for l=1:L
    mans.(['M', num2str(l)]) = spherefactory(2*l+1, 1);
end
M = productmanifold(mans);

problem = struct();
problem.M = M;

opts = struct();
opts.verbosity = 1;

[bispInds, symComp] = utilityFunc(L);

% Test results struct
results = struct();

% The test
for trial=1:trialsNo
    shc = r2c(randomNormalizedSHC(L, 0), L);
    shc(1) = 1;
    shcBisp = calculateBispectrum_M(shc, L);
    
    % Inverting the bipsectrm of shc
    t = tic();
    r = Inf;
    tries = 0;
    bestResidual = Inf;
    bestSHC = shc;
    while tries<10
        problem.cost = @(x) inversionObjectiveFuncManopt(x, shcBisp, L, bispInds, symComp);
        problem.egrad = @(x) inversionObjectiveFuncGradManopt(x, shcBisp, L, bispInds, symComp);
        
        x0 = M.rand();
%         checkgradient(problem);
        [invertedSHC, rootedResidual, output] = trustregions(problem, x0, opts);
        
        
        r = output.gradnorm;
        tries = tries + 1;
        
        if bestResidual>sqrt(rootedResidual)
            bestRediaul = sqrt(rootedResidual);
            bestSHC = SHCProdFactor(r2c(struct2vecManopt(invertedSHC, L), L), 1/sqrt(2), L);
            bestOutput = output;
            initResidual = sqrt(problem.cost(x0));
        end
        
        if sqrt(rootedResidual)<10^-8 && r<=tol
            break;
        end
    end
    invertedSHC = bestSHC;
    rootedResidual = bestRediaul;
    inversionRuntime = toc(t);
    
    % Aligning the inverted SHC and the original one
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
    results(trial).initResidual = initResidual;
    results(trial).optimalRotation = optimalRotation;
    results(trial).refinedMaxCorrelation = refinedMaxCorrelation;
    results(trial).crudeMaxCorrelation = crudeMaxCorrelation;
    results(trial).fminserachOutput = fminserachOutput;
    results(trial).inversionRuntime = inversionRuntime;
    results(trial).alignmentRuntime = alignmentRuntime;
    results(trial).tries = tries;
    results(trial).maxTries = maxTries;
    results(trial).bestRediaul = bestRediaul;
    results(trial).bestSHC = bestSHC;
    results(trial).bestOutput = bestOutput;
    
    
    disp(['iter #', num2str(trial), ' of ', num2str(trialsNo), ...
            '. Bandlimit = ', num2str(L), ...
            '. Inv runtime = ', num2str(inversionRuntime), ...
            '. Align runtime = ', num2str(alignmentRuntime), '.']);
    disp(['init res = ', num2str(initResidual), ...
        '. final res = ', num2str(rootedResidual), ...
        '. rel dist = ', num2str(relativeDistance)]);
end

Task7L6 = results;
save('Task7L6.mat', 'Task7L6');

%%
% ===================
% Ploting the results
% ===================

results = load('Task7L6.mat');
results = results.Task7L6;

color =     [0, 0.4470, 0.7410; ...
            [0.8500 0.3250 0.0980]; ...
            [0.9290 0.6940 0.1250]; ...
            [0.4940 0.1840 0.5560]; ...
            [0.4660 0.6740 0.1880]];
sz = 3;

fig = figure;
tiledlayout(3, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile;
hold;
y = [results.relativeDistance];
y = y(:);
plot(y, ...
    'Marker', 'o', ...
    'MarkerSize', sz);
hold off;
set(gca, 'yscale', 'log');

title('Bispectrum inversion of real-valued spherical function');
ylabel('Relative error');
% yticks(0:0.1:1);

% ============

nexttile;
hold;
J = 1;
y = [results.rootedResidual];
y = y(:);
y2 = [results.initResidual];
y2 = y2(:);

plot(y2, ...
    'Marker', 'o', ...
    'MarkerSize', sz);
plot(y, ...
    'Marker', 'o', ...
    'MarkerSize', sz);
hold off;
set(gca, 'yscale', 'log');
ylabel('Residual (rooted)');

legend({'Initial', 'Final'}, ...
    'Orientation', 'vertical', 'Location', 'northeast');
ylim([2.6, 3.6]);
yticks(2.6:0.2:3.6);
% ============

nexttile;

hold;
N = length(results);
S = zeros(N, 1);
for J=1:N
    S(J) = results(J).output(end).gradnorm;
end

y = S;
y = y(:);
plot(y, ...
    'Marker', 'o', ...
    'MarkerSize', sz);
hold off;
set(gca, 'yscale', 'log');

xlabel('Trial');
ylabel('First order optimality condition');
% yticks(10.^(-13:1));

save2pdf('Task7L6fig.pdf', fig);

%% Task 1
% Frequency marching in the real-valued case. Does it help?

%% 
L = 5;
td = loadtd('sf025.00339');
trialsNo = 100;

opts = optimoptions(@lsqnonlin, ...
    'SpecifyObjectiveGradient', true, ...
    'OptimalityTolerance', 10^-12, ...
    'FunctionTolerance', 10^-15, ...
    'StepTolerance', 10^-10, ...
    'Display', 'off');

[bispInds, symComp] = utilityFunc(L);
bispIndsFM = utilityFuncFreqMarch(L);

tol = 10^-8;
maxTries = 10;

% Test results struct
results = struct();

% The test
for trial=1:trialsNo
    shc = r2c(randomNormalizedSHC(L, 0), L);
    shc(1) = 1;
    shcBisp = calculateBispectrum_M(shc, L);
    
    % Inverting the bipsectrm of shc
    t = tic();
    
    invertedSHC = zeros((L+1)^2, 1);
    invertedSHC(1) = nthroot(shcBisp(1), 3);
    
    for K=1:L
        t = tic();
        r = Inf;
        tries = 0;
        rootedResidual(K) = Inf;
        bestSHC = shc;

        while tries<10
            func = @(x) inversionObjectiveFuncRealFM([invertedSHC(1:K^2); x], shcBisp, L, K, bispInds, symComp, bispIndsFM);
            x0 = randomNormalizedSHC(L, 0);
            x0 = x0(K^2+1:(K+1)^2);

            [invertedSHCKMode, rootedResidualMode, ~, ~, output] = lsqnonlin(func, x0, [], [], opts);
            
            r = output.firstorderopt;
            tries = tries + 1;

            if rootedResidual(K)>sqrt(rootedResidualMode)
                rootedResidual(K) = sqrt(rootedResidualMode);
                invertedSHC((K^2+1):(K+1)^2) = invertedSHCKMode;
                allFirstOrderOpt(K) = output.firstorderopt;
            end
            
            if sqrt(rootedResidual(K))<10^-8 && r<=tol
                break;
            end
        end
    end
    inversionRuntime = toc(t);
    invertedSHC = r2c(invertedSHC, L);
    
    % Aligning the inverted SHC and the original one
    t = tic();
    [relativeDistance, alignedSHC2, optimalRotation, refinedMaxCorrelation, crudeMaxCorrelation, fminserachOutput] = alignSphericalHarmonics(shc, invertedSHC, L, td);
    alignmentRuntime = toc(t);

    % Save results
    results(trial).bandlimit = L;
    results(trial).shc = shc;
    results(trial).shcBisp = shcBisp;
    results(trial).invertedSHC = invertedSHC;
    results(trial).rootedResidual = rootedResidual;
    results(trial).allFirstOrderOpt = allFirstOrderOpt;
    results(trial).relativeDistance = relativeDistance;
    results(trial).alignedSHC2 = alignedSHC2;
    results(trial).optimalRotation = optimalRotation;
    results(trial).refinedMaxCorrelation = refinedMaxCorrelation;
    results(trial).crudeMaxCorrelation = crudeMaxCorrelation;
    results(trial).fminserachOutput = fminserachOutput;
    results(trial).inversionRuntime = inversionRuntime;
    results(trial).alignmentRuntime = alignmentRuntime;
    results(trial).tries = tries;
    results(trial).maxTries = maxTries;
    
    
    disp(['iter #', num2str(trial), ' of ', num2str(trialsNo), ...
            '. Bandlimit = ', num2str(L), ...
            '. Inv runtime = ', num2str(inversionRuntime), ...
            '. Align runtime = ', num2str(alignmentRuntime), '.']);
    disp(['final res = ', num2str(rootedResidual), ...
        '. rel dist = ', num2str(relativeDistance)]);
end

Task1L5 = results;
save('Task1L5.mat', 'Task1L5');

%%
% ===================
% Ploting the results
% ===================

results = load('Task1L5.mat');
results = results.Task1L5;

color =     [0, 0.4470, 0.7410; ...
            [0.8500 0.3250 0.0980]; ...
            [0.9290 0.6940 0.1250]; ...
            [0.4940 0.1840 0.5560]; ...
            [0.4660 0.6740 0.1880]];
sz = 3;

fig = figure;
tiledlayout(3, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile;
hold;
y = [results.relativeDistance];
y = y(:);
plot(y, ...
    'Marker', 'o', ...
    'MarkerSize', sz);
hold off;

title('Bispectrum inversion of real-valued spherical function');
ylabel('Relative error');
ylim([0.4, 0.9]);
yticks(0.4:0.1:0.9);
% ============

nexttile;
hold;
N = length(results);
y = zeros(N, 1);
for J=1:length(results)
    y(J) = max(results(J).rootedResidual);
end

plot(y, ...
    'Marker', 'o', ...
    'MarkerSize', sz);
hold off;
set(gca, 'yscale', 'log');
ylabel('Maximal residual (rooted)');

% ============

nexttile;

hold;
N = length(results);
S = zeros(N, 1);
for J=1:N
    y(J) = max(results(J).allFirstOrderOpt);
end
y = y(:);
plot(y, ...
    'Marker', 'o', ...
    'MarkerSize', sz);
hold off;
set(gca, 'yscale', 'log');

xlabel('Trial');
ylabel('Maximal first ord. opt. cond.');

save2pdf('Task1L6fig.pdf', fig);
