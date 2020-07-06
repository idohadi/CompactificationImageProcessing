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

%% test2
fig = figure();
tiledlayout(3, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

expNo = 5;
trialsNo = 500;
MrkSz = 3;

nexttile;

hold;

S = [complexValuedSHCInverstionInitialPointTest(:, 1).relativeDistance];
plot(1:trialsNo, S, ...
    'Marker', 'o', ...
    'MarkerFaceColor', [0 0.4470 0.7410], ...
    'MarkerSize', MrkSz);

S = [complexValuedSHCInverstionInitialPointTest(:, 2).relativeDistance];
plot(1:trialsNo, S, ...
    'Marker', 'o', ...
    'MarkerFaceColor', [0.8500 0.3250 0.0980], ...
    'MarkerSize', MrkSz);

S = [complexValuedSHCInverstionInitialPointTest(:, 3).relativeDistance];
plot(1:trialsNo, S, ...
    'Marker', 'o', ...
    'MarkerFaceColor', [0.9290 0.6940 0.1250], ...
    'MarkerSize', MrkSz);

S = [complexValuedSHCInverstionInitialPointTest(:, 4).relativeDistance];
plot(1:trialsNo, S, ...
    'Marker', 'o', ...
    'MarkerFaceColor', [0.4940 0.1840 0.5560], ...
    'MarkerSize', MrkSz);

S = [complexValuedSHCInverstionInitialPointTest(:, 5).relativeDistance];
plot(1:trialsNo, S, ...
    'Marker', 'o', ...
    'MarkerFaceColor', [0.4660 0.6740 0.1880], ...
    'MarkerSize', MrkSz);

hold off;

% [0 0.4470 0.7410]
% [0.8500 0.3250 0.0980]
% [0.9290 0.6940 0.1250]
% [0.4940 0.1840 0.5560]
% [0.4660 0.6740 0.1880]

set(gca, 'yscale', 'log');

title('Bispectrum inversion of complex-valued spherical function');
xlabel('Iteration');
ylabel('Relative error');

% legend({'bandlimit = 4', 'bandlimit = 5', 'bandlimit = 6', ...
%     'bandlimit = 7', 'bandlimit = 8'}, ...
%     'Location', 'southoutside', 'Orientation', 'horizontal');
% xticks(1:trialsNo);
% yticks(10.^(-6:0));

% ============

nexttile;

hold;

S = [complexValuedSHCInverstionInitialPointTest(:, 1).rootedResidual];
plot(1:trialsNo, S, ...
    'Marker', 'o', ...
    'MarkerFaceColor', [0 0.4470 0.7410], ...
    'MarkerSize', MrkSz);

S = [complexValuedSHCInverstionInitialPointTest(:, 2).rootedResidual];
plot(1:trialsNo, S, ...
    'Marker', 'o', ...
    'MarkerFaceColor', [0.8500 0.3250 0.0980], ...
    'MarkerSize', MrkSz);

S = [complexValuedSHCInverstionInitialPointTest(:, 3).rootedResidual];
plot(1:trialsNo, S, ...
    'Marker', 'o', ...
    'MarkerFaceColor', [0.9290 0.6940 0.1250], ...
    'MarkerSize', MrkSz);

S = [complexValuedSHCInverstionInitialPointTest(:, 4).rootedResidual];
plot(1:trialsNo, S, ...
    'Marker', 'o', ...
    'MarkerFaceColor', [0.4940 0.1840 0.5560], ...
    'MarkerSize', MrkSz);

S = [complexValuedSHCInverstionInitialPointTest(:, 5).rootedResidual];
plot(1:trialsNo, S, ...
    'Marker', 'o', ...
    'MarkerFaceColor', [0.4660 0.6740 0.1880], ...
    'MarkerSize', MrkSz);

hold off;

% [0 0.4470 0.7410]
% [0.8500 0.3250 0.0980]
% [0.9290 0.6940 0.1250]
% [0.4940 0.1840 0.5560]
% [0.4660 0.6740 0.1880]

set(gca, 'yscale', 'log');

% title('Bispectrum inversion of complex-valued spherical function');
xlabel('Iteration');
ylabel('Residual (rooted)');

% legend({'bandlimit = 4', 'bandlimit = 5', 'bandlimit = 6', ...
%     'bandlimit = 7', 'bandlimit = 8'}, ...
%     'Location', 'southoutside', 'Orientation', 'horizontal');
% xticks(1:trialsNo);

% ============

nexttile;


hold;

S = [complexValuedSHCInverstionInitialPointTest(:, 1).output];
S = [S.firstorderopt];
plot(1:trialsNo, S, ...
    'Marker', 'o', ...
    'MarkerFaceColor', [0 0.4470 0.7410], ...
    'MarkerSize', MrkSz);

S = [complexValuedSHCInverstionInitialPointTest(:, 2).output];
S = [S.firstorderopt];
plot(1:trialsNo, S, ...
    'Marker', 'o', ...
    'MarkerFaceColor', [0.8500 0.3250 0.0980], ...
    'MarkerSize', MrkSz);

S = [complexValuedSHCInverstionInitialPointTest(:, 3).output];
S = [S.firstorderopt];
plot(1:trialsNo, S, ...
    'Marker', 'o', ...
    'MarkerFaceColor', [0.9290 0.6940 0.1250], ...
    'MarkerSize', MrkSz);

S = [complexValuedSHCInverstionInitialPointTest(:, 4).output];
S = [S.firstorderopt];
plot(1:trialsNo, S, ...
    'Marker', 'o', ...
    'MarkerFaceColor', [0.4940 0.1840 0.5560], ...
    'MarkerSize', MrkSz);

S = [complexValuedSHCInverstionInitialPointTest(:, 5).output];
S = [S.firstorderopt];
plot(1:trialsNo, S, ...
    'Marker', 'o', ...
    'MarkerFaceColor', [0.4660 0.6740 0.1880], ...
    'MarkerSize', MrkSz);

hold off;

% [0 0.4470 0.7410]
% [0.8500 0.3250 0.0980]
% [0.9290 0.6940 0.1250]
% [0.4940 0.1840 0.5560]
% [0.4660 0.6740 0.1880]

set(gca, 'yscale', 'log');

% title('Bispectrum inversion of complex-valued spherical function');
xlabel('Iteration');
ylabel('First order optimality condition');

legend({'experiment 1', 'experiment 2', 'experiment 3', ...
    'experiment 4', 'experiment 5'}, ...
    'Location', 'southoutside', 'Orientation', 'horizontal');
% xticks(1:trialsNo);
