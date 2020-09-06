%% Parameters
Ls = [5, 6, 8, 10, 12];
color =     [0, 0.4470, 0.7410; ...
            [0.8500 0.3250 0.0980]; ...
            [0.9290 0.6940 0.1250]; ...
            [0.4940 0.1840 0.5560]; ...
            [0.4660 0.6740 0.1880]];
trialsNo = 20;

%% Test1: complex-valued function inversion

results = load('complexValuedSHCInversionTestResultsNewMEXFunc.mat');
results = results.complexValuedSHCInversionTestResultsNewMEXFunc;

fig = figure();
tiledlayout(3, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile;
hold;
for J=1:length(Ls)
    plot(1:trialsNo, [results(:, J).relativeDistance], ...
        'Marker', 'o', ...
        'MarkerFaceColor', color(J, :));
end
hold off;
set(gca, 'yscale', 'log');

title('Bispectrum inversion of complex-valued spherical function');
xlabel('Iteration');
ylabel('Relative error');
xticks(1:trialsNo);

% ============

nexttile;
hold;
for J=1:length(Ls)
    plot(1:trialsNo, [results(:, J).rootedResidual], ...
        'Marker', 'o', ...
        'MarkerFaceColor', color(J, :));
end
hold off;
set(gca, 'yscale', 'log');

xlabel('Iteration');
ylabel('Residual (rooted)');
xticks(1:trialsNo);

% ============

nexttile;
S = [results.output];
S = reshape(S, [trialsNo, length(Ls)]);
hold;
for J=1:length(Ls)
    plot(1:trialsNo, [S(:, J).firstorderopt], ...
        'Marker', 'o', ...
        'MarkerFaceColor', color(J, :));
end
hold off;

set(gca, 'yscale', 'log');

xlabel('Iteration');
ylabel('First order optimality condition');

legend({'bandlimit = 4', 'bandlimit = 6', 'bandlimit = 8', ...
    'bandlimit = 10', 'bandlimit = 12'}, ...
    'Location', 'southoutside', 'Orientation', 'horizontal');
xticks(1:trialsNo);


%% Test1 FOLLOWUP: complex-valued function inversion
results = load('complexValuedSHCInversionTestResultsNewMEXFuncFollowup.mat');
results = results.complexValuedSHCInversionTestResultsNewMEXFuncFollowup;

fig = figure();
tiledlayout(3, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile;
hold;
for J=1:length(results)
    plot(1:length([results{J}.relativeDistance]), [results{J}.relativeDistance], ...
        'Marker', 'o', ...
        'MarkerFaceColor', color(J, :));
end
hold off;
set(gca, 'yscale', 'log');

title('Bispectrum inversion of complex-valued spherical function');
xlabel('Iteration');
ylabel('Relative error');
xticks(1:trialsNo);

% ============

nexttile;
hold;
for J=1:length(results)
    plot(1:length(results{J}), [results{J}.rootedResidual], ...
        'Marker', 'o', ...
        'MarkerFaceColor', color(J, :));
end
hold off;
set(gca, 'yscale', 'log');

xlabel('Iteration');
ylabel('Residual (rooted)');
xticks(1:trialsNo);

% ============

nexttile;
hold;
for J=1:length(Ls)
    S = [results{J}.output];
    plot(1:length(S), [S.firstorderopt], ...
        'Marker', 'o', ...
        'MarkerFaceColor', color(J, :));
end
hold off;

set(gca, 'yscale', 'log');

xlabel('Iteration');
ylabel('First order optimality condition');

legend({'bandlimit = 4', 'bandlimit = 6', 'bandlimit = 8', ...
    'bandlimit = 10', 'bandlimit = 12'}, ...
    'Location', 'southoutside', 'Orientation', 'horizontal');
xticks(1:trialsNo);


%% Test 2: real-valued function inversion
Ls = [5, 6, 8, 10, 12];
color =     [0, 0.4470, 0.7410; ...
            [0.8500 0.3250 0.0980]; ...
            [0.9290 0.6940 0.1250]; ...
            [0.4940 0.1840 0.5560]; ...
            [0.4660 0.6740 0.1880]];
trialsNo = 20;

results = load('realValuedSHCInversionTestResultsNewMEXFunc.mat');
results = results.realValuedSHCInversionTestResultsNewMEXFunc;

fig = figure();
tiledlayout(3, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile;
hold;
for J=1:size(results, 2)
    plot(1:trialsNo, [results(:, J).relativeDistance], ...
        'Marker', 'o', ...
        'MarkerFaceColor', color(J, :));
end
hold off;

title('Bispectrum inversion of real-valued spherical function');
xlabel('Iteration');
ylabel('Relative error');
xticks(1:trialsNo);
ylim([0.2, 1]);

% ============

nexttile;
hold;
for J=1:size(results, 2)
    plot(1:trialsNo, [results(:, J).rootedResidual], ...
        'Marker', 'o', ...
        'MarkerFaceColor', color(J, :));
end
hold off;
set(gca, 'yscale', 'log');

xlabel('Iteration');
ylabel('Residual (rooted)');
xticks(1:trialsNo);
yticks(10.^(-10:-1));

% ============

nexttile;
hold;
for J=1:size(results, 2)
    S = [results(:, J).output];
    plot(1:trialsNo, [S.firstorderopt], ...
        'Marker', 'o', ...
        'MarkerFaceColor', color(J, :));
end
hold off;

set(gca, 'yscale', 'log');

xlabel('Iteration');
ylabel('First order optimality condition');
yticks(10.^(-10:-1));

legcell = {};
for J=1:size(results, 2)
    legcell{J} = ['bandlimit = ', num2str(results(1, J).bandlimit)];
end
legend(legcell, ...
    'Location', 'southoutside', 'Orientation', 'horizontal');
xticks(1:trialsNo);