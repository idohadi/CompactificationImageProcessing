%% Real-valued bispectrum inversion for large bandlimit
%% 
L = 12;
td = loadtd('sf025.00339');
trialsNo = 30;

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
    func = @(x) inversionObjectiveFunc(cfy(x), shcBisp, L);
    % Generating initial point
    x0 = r2c(randomNormalizedSHC(L, 0), L);
    % Inverting the bipsectrm of shc
    [invertedSHC, rootedResidual, ~, ~, output] = lsqnonlin(func, rlfy(x0), [], [], opts);
    invertedSHC = cfy(invertedSHC);
    rootedResidual = sqrt(rootedResidual);
    inversionRuntime = toc(t);
    
    initDistance = norm(x0-shc, 2);
    initRelativeDistance = initDistance/norm(shc, 2);

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
    results(trial).x0 = x0;        
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
    results(trial).initRelativeDistance = initRelativeDistance;

    disp(['trial #', num2str(trial), ' of ', num2str(trialsNo), ...
        '. Bandlimit = ', num2str(L), ...
        '. Inv runtime = ', num2str(inversionRuntime), ...
        '. Align runtime = ', num2str(alignmentRuntime), '.']);
    disp(['. init dist = ', num2str(initDistance), ...
        '. init rel dist = ', num2str(initRelativeDistance), ...
        '. final res = ', num2str(rootedResidual), ...
        '. rel dist = ', num2str(relativeDistance), ...
        '. first ord. opt. = ', num2str(output.firstorderopt)]);
end

LargeBandlimitL12 = results;
save('LargeBandlimitL12.mat', 'LargeBandlimitL12');

%%
% ===================
% Ploting the results
% ===================

results = load('LargeBandlimitL12.mat');
results = results.LargeBandlimitL12;

color =     [0, 0.4470, 0.7410; ...
            [0.8500 0.3250 0.0980]; ...
            [0.9290 0.6940 0.1250]; ...
            [0.4940 0.1840 0.5560]; ...
            [0.4660 0.6740 0.1880]];
sz = 3;

fig = figure;
tiledlayout(4, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

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
ylabel('Relative error');
yticks(0:0.05:1.5);

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

yticks(0:0.05:1.5);
ylabel('Residual (rooted)');

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

ylabel('First ord. opt. cond.');
yticks(10.^(-13:1));


nexttile;

hold;
x = [results.initDistance];
x = x(:);
[x, inds] = sortrows(x);
y1 = [results.inversionRuntime];
y1 = y1(:);
y1 = y1(inds);
y2 = [results.alignmentRuntime];
y2 = y2(:);
y2 = y2(inds);
plot(x, y1, ...
    'Marker', 'o', ...
    'MarkerSize', sz);
plot(x, y2, ...
    'Marker', 'o', ...
    'MarkerSize', sz);
hold off;
legend({'Inversion', 'Alignment'}, ...
    'Location', 'northeast', ...
    'Orientation', 'vertical');
ylabel('Runtime (sec)');
xlim([min(x), max(x)]);
xlabel('Distance between initial guess and ground truth');

save2pdf('LargeBandlimitL12fig.pdf', fig);

%% Task 2 for large bandlimit
% In the real-valued case, what is the effect of the distance between the 
% initial point and the ground truth?
%% 
L = 12;
td = loadtd('sf025.00339');
trialsNo = 30;
pow_lb = -5;
pow_ub = 2;

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
    func = @(x) inversionObjectiveFunc(cfy(x), shcBisp, L);
    % Generating initial point
    direction = r2c(randomNormalizedSHC(L, 0), L);
    direction = direction/norm(direction, 2);
    digit = randi(9);
    pow = randi(pow_ub - pow_lb + 1) + pow_lb - 1;
    initDistance = digit*10^pow;
    x0 = initDistance*direction + shc;
    % Inverting the bipsectrm of shc
    [invertedSHC, rootedResidual, ~, ~, output] = lsqnonlin(func, rlfy(x0), [], [], opts);
    invertedSHC = cfy(invertedSHC);
    rootedResidual = sqrt(rootedResidual);
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
    results(trial).x0 = x0;        
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
    disp(['. init dist = ', num2str(initDistance), ...
        '. final res = ', num2str(rootedResidual), ...
        '. rel dist = ', num2str(relativeDistance), ...
        '. first ord. opt. = ', num2str(output.firstorderopt)]);
end

Task2L12 = results;
save('Task2L12.mat', 'Task2L12');

%%
% ===================
% Ploting the results
% ===================
% TODO: write this code

results = load('Task2L12.mat');
results = results.Task2L12;
trialsNo = length(results);

color =     [0, 0.4470, 0.7410; ...
            [0.8500 0.3250 0.0980]; ...
            [0.9290 0.6940 0.1250]; ...
            [0.4940 0.1840 0.5560]; ...
            [0.4660 0.6740 0.1880]];
sz = 3;

fig = figure;
tiledlayout(4, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

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
set(gca, 'yscale', 'log');

title('Bispectrum inversion of real-valued spherical function');
ylabel('Relative error');
yticks(10.^(-13:1));

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

set(gca, 'yscale', 'log');
set(gca, 'xscale', 'log');

yticks(10.^(-13:1));
ylabel('Residual (rooted)');

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

ylabel('First ord. opt. cond.');
yticks(10.^(-13:1));


nexttile;

hold;
x = [results.initDistance];
x = x(:);
[x, inds] = sortrows(x);
y1 = [results.inversionRuntime];
y1 = y1(:);
y1 = y1(inds);
y2 = [results.alignmentRuntime];
y2 = y2(:);
y2 = y2(inds);
plot(x, y1, ...
    'Marker', 'o', ...
    'MarkerSize', sz);
plot(x, y2, ...
    'Marker', 'o', ...
    'MarkerSize', sz);
hold off;
legend({'Inversion', 'Alignment'}, ...
    'Location', 'southeast', ...
    'Orientation', 'vertical');
ylabel('Runtime (sec)');
set(gca, 'xscale', 'log');
set(gca, 'yscale', 'log');
% ylim([10^0, 10^2]);
yticks(10.^(-13:5));

xlabel('Distance between initial guess and ground truth');


save2pdf('Task2L12fig.pdf', fig);


%% Task 2 for several large bandlimit
% In the real-valued case, what is the effect of the distance between the 
% initial point and the ground truth?
%% 
Ls = [8, 10, 12, 14, 16];
tds = ['sf018.00182'; 'sf022.00266'; 'sf026.00366'; 'sf030.00482'; 'sf034.00614'];

trialsNo = 30;
% Lower and upper relative error in percents
rel_err_lb = 5;
rel_err_ub = 100;

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
for J=1:length(Ls)
    L = Ls(J);
    td = loadtd(tds(J, :));

    for trial=1:trialsNo
        shc = r2c(randomNormalizedSHC(L, 0), L);
        shcBisp = calculateBispectrum_M(shc, L);

        t = tic();
        func = @(x) inversionObjectiveFunc(cfy(x), shcBisp, L);
        % Generating initial point
        direction = r2c(randomNormalizedSHC(L, 0), L);
        direction = direction/norm(direction, 2);
        initialRelativeErrorPercent = (rel_err_ub-rel_err_lb)*rand() + rel_err_lb;
        initialRelativeError = initialRelativeErrorPercent/100;
        initDistance = sqrt(L+1)*initialRelativeError;
        x0 = initDistance*direction + shc;
        % Inverting the bipsectrm of shc
        [invertedSHC, rootedResidual, ~, ~, output] = lsqnonlin(func, rlfy(x0), [], [], opts);
        invertedSHC = cfy(invertedSHC);
        rootedResidual = sqrt(rootedResidual);
        inversionRuntime = toc(t);


        % Aligning the inverted SHC and the original one
        t = tic();
        [relativeDistance, alignedSHC2, optimalRotation, refinedMaxCorrelation, crudeMaxCorrelation, fminserachOutput] = alignSphericalHarmonics(shc, invertedSHC, L, td);
        alignmentRuntime = toc(t);

        % Save results
        results(trial, J).bandlimit = L;
        results(trial, J).shc = shc;
        results(trial, J).shcBisp = shcBisp;
        
        results(trial, J).direction = direction;
        results(trial, J).initDistance = initDistance;
        results(trial, J).initialRelativeErrorPercent = initialRelativeErrorPercent;
        results(trial, J).initialRelativeError = initialRelativeError;
        results(trial, J).x0 = x0;
                
        results(trial, J).invertedSHC = invertedSHC;
        results(trial, J).rootedResidual = rootedResidual;
        results(trial, J).output = output;
        results(trial, J).relativeDistance = relativeDistance;
        results(trial, J).alignedSHC2 = alignedSHC2;
        results(trial, J).optimalRotation = optimalRotation;
        results(trial, J).refinedMaxCorrelation = refinedMaxCorrelation;
        results(trial, J).crudeMaxCorrelation = crudeMaxCorrelation;
        results(trial, J).fminserachOutput = fminserachOutput;
        results(trial, J).inversionRuntime = inversionRuntime;
        results(trial, J).alignmentRuntime = alignmentRuntime;

        disp(['iter #', num2str(trial), ' of ', num2str(trialsNo), ...
            '. Bandlimit = ', num2str(L), ...
            '. Inv runtime = ', num2str(inversionRuntime), ...
            '. Align runtime = ', num2str(alignmentRuntime), '.']);
        disp(['. rel init dist = ', num2str(initialRelativeError), ...
            '. final res = ', num2str(rootedResidual), ...
            '. rel dist = ', num2str(relativeDistance), ...
            '. first ord. opt. = ', num2str(output.firstorderopt)]);
    end
end

Task2SeveralLargeLAllExpanded = results;
save('Task2SeveralLargeLAllExpanded.mat', 'Task2SeveralLargeLAllExpanded');

%%
% ===================
% Ploting the results
% ===================

results = load('Task2SeveralLargeLAllExpanded.mat');
results = results.Task2SeveralLargeLAllExpanded;

sz = 3;

fig = figure;
tiledlayout(5, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile;
hold;
for J=1:size(results, 2)
    x = [results(:, J).initialRelativeError];
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
set(gca, 'yscale', 'log');

title('Bispectrum inversion of real-valued spherical function');
ylabel('Relative error');
yticks(10.^(-13:3));

% ============

nexttile;
hold;
for J=1:size(results, 2)
    x = [results(:, J).initialRelativeError];
    x = x(:);
    [x, inds] = sortrows(x);
    y = [results(:, J).rootedResidual];
    y = y(:);
    y = y(inds);
    plot(x, y, ...
        'Marker', 'o', ...
        'MarkerSize', sz);
end
hold off;
set(gca, 'yscale', 'log');

ylabel('Residual (rooted)');
yticks(10.^(-13:3));

% ============

nexttile;

hold;
for J=1:size(results, 2)
    S = [results(:, J).output];

    x = [results(:, J).initialRelativeError];
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

ylabel('First ord. opt. cond.');
yticks(10.^(-13:3));


nexttile;

hold;
for J=1:size(results, 2)
    x = [results(:, J).initialRelativeError];
    x = x(:);
    [x, inds] = sortrows(x);
    y = [results(:, J).inversionRuntime];
    y = y(:);
    y = y(inds);
    plot(x, y, ...
        'Marker', 'o', ...
        'MarkerSize', sz);
end
hold off;
set(gca, 'yscale', 'log');

ylabel('Inversion untime (sec)');
yticks(10.^(-13:3));

nexttile;

hold;
for J=1:size(results, 2)
    x = [results(:, J).initialRelativeError];
    x = x(:);
    [x, inds] = sortrows(x);
    y = [results(:, J).alignmentRuntime];
    y = y(:);
    y = y(inds);
    plot(x, y, ...
        'Marker', 'o', ...
        'MarkerSize', sz);
end
hold off;
set(gca, 'yscale', 'log');

ylabel('Alignment runtime (sec)');

legcell = cell(1, size(results, 2));
for J=1:size(results, 2)
    legcell{J} = ['bandlimit = ', num2str(results(1, J).bandlimit)];
end
legend(legcell, ...
    'Orientation', 'horizontal', ...
    'Location', 'southoutside');

xlabel('Initial error relative to ground truth');


save2pdf('Task2SeveralLargeLAllfig.pdf', fig);

%% Bispectrum inversion stability
% In the real-valued case, what is the effect of the distance between the 
% initial point and the ground truth?
%% 
Ls = [8, 10, 12, 14, 16];
tds = ['sf018.00182'; 'sf022.00266'; 'sf026.00366'; 'sf030.00482'; 'sf034.00614'];

trialsNo = 50;
% Lower and upper relative error in percents
rel_err_lb = 10^-2;
rel_err_ub = 15;

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
for J=1:length(Ls)
    L = Ls(J);
    td = loadtd(tds(J, :));
    
    [bispInds, ~] = utilityFunc(L);
    
    for trial=1:trialsNo
        shc = r2c(randomNormalizedSHC(L, 0), L);
        shcBisp = calculateBispectrum_M(shc, L);
        % Perturb the bispectrum
        shcBispPerturbationDirection = randn(length(shcBisp), 1);
        shcBispPerturbationDirection(~bispInds) = 0;
        shcBispPerturbationDirection ...
            = shcBispPerturbationDirection/norm(shcBispPerturbationDirection, 2);
        initialBispRelativeErrorPercent = (rel_err_ub-rel_err_lb)*rand() + rel_err_lb;
        initialBispRelativeError = initialBispRelativeErrorPercent/100;
        initBispDistance = norm(shcBisp, 2)*initialBispRelativeError;
        shcBispPerturbed = shcBisp + initBispDistance*shcBispPerturbationDirection;
        
        t = tic();
        func = @(x) inversionObjectiveFunc(cfy(x), shcBispPerturbed, L);
        % Generating initial point
        direction = r2c(randomNormalizedSHC(L, 0), L);
        direction = direction/norm(direction, 2);
        initialRelativeErrorPercent = 60;
        initialRelativeError = initialRelativeErrorPercent/100;
        initDistance = sqrt(L+1)*initialRelativeError;
        x0 = initDistance*direction + shc;
        % Inverting the bipsectrm of shc
        [invertedSHC, rootedResidual, ~, ~, output] = lsqnonlin(func, rlfy(x0), [], [], opts);
        invertedSHC = cfy(invertedSHC);
        rootedResidual = sqrt(rootedResidual);
        inversionRuntime = toc(t);


        % Aligning the inverted SHC and the original one
        t = tic();
        [relativeDistance, alignedSHC2, optimalRotation, refinedMaxCorrelation, crudeMaxCorrelation, fminserachOutput] = alignSphericalHarmonics(shc, invertedSHC, L, td);
        alignmentRuntime = toc(t);

        % Save results
        results(trial, J).bandlimit = L;
        results(trial, J).shc = shc;
        results(trial, J).shcBisp = shcBisp;
        
        results(trial, J).shcBispPerturbationDirection = shcBispPerturbationDirection;
        results(trial, J).initialBispRelativeErrorPercent = initialBispRelativeErrorPercent;
        results(trial, J).initialBispRelativeError = initialBispRelativeError;
        results(trial, J).initBispDistance = initBispDistance;
        results(trial, J).shcBispPerturbed = shcBispPerturbed;
        
        results(trial, J).direction = direction;
        results(trial, J).initDistance = initDistance;
        results(trial, J).initialRelativeErrorPercent = initialRelativeErrorPercent;
        results(trial, J).initialRelativeError = initialRelativeError;
        results(trial, J).x0 = x0;
                
        results(trial, J).invertedSHC = invertedSHC;
        results(trial, J).rootedResidual = rootedResidual;
        results(trial, J).output = output;
        results(trial, J).relativeDistance = relativeDistance;
        results(trial, J).alignedSHC2 = alignedSHC2;
        results(trial, J).optimalRotation = optimalRotation;
        results(trial, J).refinedMaxCorrelation = refinedMaxCorrelation;
        results(trial, J).crudeMaxCorrelation = crudeMaxCorrelation;
        results(trial, J).fminserachOutput = fminserachOutput;
        results(trial, J).inversionRuntime = inversionRuntime;
        results(trial, J).alignmentRuntime = alignmentRuntime;

        disp(['iter #', num2str(trial), ' of ', num2str(trialsNo), ...
            '. Bandlimit = ', num2str(L), ...
            '. Inv runtime = ', num2str(inversionRuntime), ...
            '. Align runtime = ', num2str(alignmentRuntime), '.']);
        disp(['. rel init dist = ', num2str(initialRelativeError), ...
            '. rel bisp pert size = ', num2str(initialBispRelativeError), ...
            '. final res = ', num2str(rootedResidual), ...
            '. rel dist = ', num2str(relativeDistance), ...
            '. first ord. opt. = ', num2str(output.firstorderopt)]);
    end
end

RealValuedBispInvStability = results;
save('RealValuedBispInvStability.mat', 'RealValuedBispInvStability');

%%
% ===================
% Ploting the results
% ===================

results = load('RealValuedBispInvStability.mat');
results = results.RealValuedBispInvStability;

sz = 3;

fig = figure;
tiledlayout(3, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile;
hold;
for J=1:size(results, 2)
    x = [results(:, J).initialBispRelativeError];
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

title('Bispectrum inversion of real-valued spherical function');
ylabel('Relative error');
yticks(0:0.1:1);

% ============

nexttile;
hold;
for J=1:size(results, 2)
    x = [results(:, J).initialBispRelativeError];
    x = x(:);
    [x, inds] = sortrows(x);
    y = [results(:, J).rootedResidual] - [results(:, J).initBispDistance];
    y = y(:);
    y = y(inds);
    plot(x, y, ...
        'Marker', 'o', ...
        'MarkerSize', sz);
end
hold off;

ylabel('Residual, final-truth');

% ============

nexttile;

hold;
for J=1:size(results, 2)
    S = [results(:, J).output];

    x = [results(:, J).initialBispRelativeError];
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

ylabel('First ord. opt. cond.');
yticks(10.^(-13:3));

legcell = cell(1, size(results, 2));
for J=1:size(results, 2)
    legcell{J} = ['bandlimit = ', num2str(results(1, J).bandlimit)];
end
legend(legcell, ...
    'Orientation', 'horizontal', ...
    'Location', 'southoutside');

xlabel('Relative bispectrum perturbation size');


save2pdf('RealValuedBispInvStability.pdf', fig);


disp('x is relative SHC distance. y is bispectrum relative perturbation.');
for J=1:size(results, 2)
    x = [results(:, J).initialBispRelativeError];
    x = x(:);
    
    y = [results(:, J).relativeDistance];
    y = y(:);
    
    p = polyfit(x, y, 1);
    M = max(y./x, [], 'all');
    if p(2)>=0
        disp(['Bandlimit = ', num2str(results(1, J).bandlimit), ...
            '. Fitting polyonmial is ', num2str(p(1)), '*x+', num2str(p(2)), '. ', ...
            'Lip. const ~ ', num2str(M), '.']);
    else
        disp(['Bandlimit = ', num2str(results(1, J).bandlimit), ...
            '. Fitting polyonmial is ', num2str(p(1)), '*x', num2str(p(2)), '. ', ...
            'Lip. const ~ ', num2str(M), '.']);
    end
end
