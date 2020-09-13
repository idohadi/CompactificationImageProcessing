%% Init
L = 8;
CGs = build_CGs_vec(L);

%% "Truth"
[f, f_bisp] = randComplexSHC(L, CGs);

%%
f2 = f;
f2(1) = 0;
f_bisp2 = bisp_vec(f2, L, CGs);

%% Parameters
% N = [100; 500; 1000; 1500; 2000; 3000];
% N = [100; 500; 1000; 1500; 2000; 3000; 4000; 5000];
N = 500:300:5000;
N = N(:);
sigma = 10;

exp_num = 20;

DenoVec = buildDenoVec(L);

%% Experimental results
f00_test = zeros(length(N), 1);
bisp_test = zeros(length(N), 1);

%% Test
t1 = tic();
for n=1:length(N)
    disp(strcat(['Exeriment number ', num2str(n), ' of ', num2str(length(N))]));
    f00_results = zeros(exp_num, 1);
    bisp_results = zeros(exp_num, 1);
    for k=1:exp_num
        shc_samples = noisySamples(f, L, N(n), sigma);
        
        f00_est = estimatef00(shc_samples);
        f00_results(k) = relmaxdist(f00_est, f(1));
        
        bisp_samples = samplesBispectrum(shc_samples, L, CGs, f00_est);
%         bisp_est = estimateBispectrum(bisp_samples, L, f00_results(k), sigma, DenoVec);
        bisp_est = estimateBispectrum(bisp_samples, L, 0, sigma, DenoVec);
        bisp_results(k) = relmaxdist(bisp_est, f_bisp2);
    end
    
    f00_test(n) = mean(f00_results);
    bisp_test(n) = mean(bisp_results);
end
t2 = toc(t1);

%%
w = 650;
h = 2*250;
factor = 1.5;
figure('Renderer', 'painters', 'Position', [10, 10, factor*w, factor*h]);


subplot(2, 2, 1);
scatter(N, f00_test); 
set(gca, 'xscale', 'log');
set(gca, 'yscale', 'log');
P_f00 = polyfit(log(N), log(f00_test), 1);
title('$\hat{f}_{0,0}$ estimator ($\tilde{f}_{0,0}$)', ...
    'FontSize', 18, ...
    'interpreter', 'latex');
xlabel('Sample size, N', ...
    'FontSize', 13, ...
    'interpreter', 'latex');
ylabel('$ \left| \tilde{f}_{0,0} - \hat{f}_{0,0} \right| / \left| \hat{f}_{0,0} \right| $', ...
    'FontSize', 15, ...
    'interpreter', 'latex');
xlim([N(1)-10, N(end)+1]);

subplot(2, 2, 2);
scatter(N, bisp_test); 
set(gca, 'xscale', 'log');
set(gca, 'yscale', 'log');
P_bisp = polyfit(log(N), log(bisp_test), 1);
title('Bispectrum estimator ($\tilde{b}_{\ell_1, \ell_2, \ell} (\textbf{f} - \tilde{f}_{0,0})$)', ...
    'FontSize', 18, ...
    'interpreter', 'latex');
xlabel('Sample size, N', ...
    'FontSize', 13, ...
    'interpreter', 'latex');
ylabel('$ \frac{\max\limits_{\ell_1,\ell_2, \ell} \left| \tilde{b}_{\ell_1, \ell_2,\ell} \left( \textbf{f}-\hat{f}_{0,0} \right) - {b}_{\ell_1, \ell_2,\ell} \left( \textbf{f}-\hat{f}_{0,0} \right) \right|}{ \max\limits_{\ell_1,\ell_2, \ell}  \left| {b}_{\ell_1, \ell_2,\ell} \left( \textbf{f}-\hat{f}_{0,0} \right) \right|}  $', ...
    'FontSize', 15, ...
    'interpreter', 'latex');
xlim([N(1)-10, N(end)+1]);

subplot(2, 2, 3);
allOneString = sprintf('%.0f,' , N);
allOneString = allOneString(1:end-1);% strip final comma

text(0.1, 0.7, {'Signal type: complex', ...
    ['$\sigma = ', num2str(sigma), '$'], ...
    ['$N=', allOneString, '$'], ...
    ['Trials per sample size: ', num2str(exp_num)], ...
    ['$\hat{f}_{0,0}$ estimation slope (log-log): ', num2str(P_f00(1))], ...
    ['Bispectrum estimation slope (log-log): ', num2str(P_bisp(1))]}, ...
    'FontSize', 13, ...
    'interpreter', 'latex');
axis off;