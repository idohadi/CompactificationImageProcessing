%% Test1
fig = figure();
tiledlayout(3, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

trialsNo = 20;


nexttile;

hold;
plot(1:trialsNo, [complexValuedSHCInversionTestResults.relativeDistance], ...
    'Marker', 'o', ...
    'MarkerFaceColor', [0 0.4470 0.7410]);
plot(1:trialsNo, [complexValuedSHCInversionTestResults2.relativeDistance], ...
    'Marker', 'o', ...
    'MarkerFaceColor', [0.8500 0.3250 0.0980]);
plot(1:trialsNo, [complexValuedSHCInversionTestResults6.relativeDistance], ...
    'Marker', 'o', ...
    'MarkerFaceColor', [0.9290 0.6940 0.1250]);
plot(1:trialsNo, [complexValuedSHCInversionTestResults7.relativeDistance], ...
    'Marker', 'o', ...
    'MarkerFaceColor', [0.4940 0.1840 0.5560]);
plot(1:trialsNo, [complexValuedSHCInversionTestResults8.relativeDistance], ...
    'Marker', 'o', ...
    'MarkerFaceColor', [0.4660 0.6740 0.1880]);
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
xticks(1:trialsNo);
% yticks(10.^(-6:0));

% ============

nexttile;

hold;
plot(1:trialsNo, [complexValuedSHCInversionTestResults.rootedResidual], ...
    'Marker', 'o', ...
    'MarkerFaceColor', [0 0.4470 0.7410]);
plot(1:trialsNo, [complexValuedSHCInversionTestResults2.rootedResidual], ...
    'Marker', 'o', ...
    'MarkerFaceColor', [0.8500 0.3250 0.0980]);
plot(1:trialsNo, [complexValuedSHCInversionTestResults6.rootedResidual], ...
    'Marker', 'o', ...
    'MarkerFaceColor', [0.9290 0.6940 0.1250]);
plot(1:trialsNo, [complexValuedSHCInversionTestResults7.rootedResidual], ...
    'Marker', 'o', ...
    'MarkerFaceColor', [0.4940 0.1840 0.5560]);
plot(1:trialsNo, [complexValuedSHCInversionTestResults8.rootedResidual], ...
    'Marker', 'o', ...
    'MarkerFaceColor', [0.4660 0.6740 0.1880]);
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
xticks(1:trialsNo);

% ============

nexttile;

hold;

S = [complexValuedSHCInversionTestResults.output];
plot(1:trialsNo, [S.firstorderopt], ...
    'Marker', 'o', ...
    'MarkerFaceColor', [0 0.4470 0.7410]);

S = [complexValuedSHCInversionTestResults2.output];
plot(1:trialsNo, [S.firstorderopt], ...
    'Marker', 'o', ...
    'MarkerFaceColor', [0.8500 0.3250 0.0980]);

S = [complexValuedSHCInversionTestResults6.output];
plot(1:trialsNo, [S.firstorderopt], ...
    'Marker', 'o', ...
    'MarkerFaceColor', [0.9290 0.6940 0.1250]);

S = [complexValuedSHCInversionTestResults7.output];
plot(1:trialsNo, [S.firstorderopt], ...
    'Marker', 'o', ...
    'MarkerFaceColor', [0.4940 0.1840 0.5560]);

S = [complexValuedSHCInversionTestResults8.output];
plot(1:trialsNo, [S.firstorderopt], ...
    'Marker', 'o', ...
    'MarkerFaceColor', [0.4660 0.6740 0.1880]);

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

legend({'bandlimit = 4', 'bandlimit = 5', 'bandlimit = 6', ...
    'bandlimit = 7', 'bandlimit = 8'}, ...
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
