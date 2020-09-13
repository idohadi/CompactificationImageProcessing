%% Init
L = 8;
CGs = build_CGs_vec(L);

%% "Truth"
[f, ~] = randComplexSHC(L, CGs);

% The bispectrum of f - f_{0,0}
f2 = f;
f2(1) = 0;
f_bisp = bisp_vec(f2, L, CGs);

%% Test parameters
% Sample sizes
% N = 100:100:500;
N = [500, 1000, 1500, 2500, 3000, 3500, 4000, 5000, 7500, 10^4, 15000, 20000];
N = N(:);

% Noise level
% sigma = [0.5; 1];
sigma = [0.5, 1, 2, 4, 8, 10];
sigma = sigma(:);

% Number of trials per experiment
trial_num = 20;

% Whether to add random rotations to samples
rotate = true;

%% Experiments
Exps_results = struct();

t1 = tic();
for s=1:length(sigma)
    disp(strcat(['Sigma = ', num2str(sigma(s))]));
    
    for n=1:length(N)
        disp(strcat(['Exeriment number ', num2str(n), ' of ', num2str(length(N)), '. N = ', num2str(N(n))]));
        
        f00_results = zeros(trial_num, 1);
        bisp_results = zeros(trial_num, 1);
        
        for k=1:trial_num
            shc_samples = noisySamples(f, L, N(n), sigma(s), rotate);
            
            f00_est = estimatef00(shc_samples);
            f00_results(k) = relmaxdist(f00_est, f(1));
            
            bisp_est = estimateBispectrumInSitu(shc_samples, L, CGs, f00_est);
            bisp_results(k) = relmaxdist(bisp_est, f_bisp);
        end
        
        Exps_results(s, n).sigma = sigma(s);
        Exps_results(s, n).N = N(n);
        Exps_results(s, n).trial_num = trial_num;
        
        Exps_results(s, n).f00_results = f00_results;
        Exps_results(s, n).bisp_results = bisp_results;
        
        Exps_results(s, n).f00_mean = mean(f00_results);
        Exps_results(s, n).bisp_mean = mean(bisp_results);
    end
end
t2 = toc(t1);


%% Ploting things
figure;

for J=1:length(sigma)
    
    subplot(3, 2, J);

    hold;

    f00_test = [Exps_results([Exps_results.sigma] == sigma(J)).f00_mean];
    bisp_test = [Exps_results([Exps_results.sigma] == sigma(J)).bisp_mean];
    scatter(N, f00_test, 12, 'filled', 'red'); 
    scatter(N, bisp_test, 12, 'filled', 'blue'); 

    P_f00 = polyfit(log(N), log(f00_test).', 1);
    P_bisp = polyfit(log(N), log(bisp_test).', 1);

    xlabel('$N$, sample size', ...
        'FontSize', 14, ...
        'interpreter', 'latex');
    ylabel({'max-norm relative', 'to ground truth'}, ...
        'FontSize', 14, ...
        'interpreter', 'latex');
    caption = strcat(['$\sigma^2 = ', num2str(sigma(J)), '$']);
    title({caption, ['$f_{0,0}$ slope: ', num2str(P_f00(1))], ['$b$ slope: ', num2str(P_bisp(1))]}, ...
        'Fontsize', 13, ...
        'interpreter', 'latex');


    xl = xlim;
    yl = ylim;
    xt = 0.06 * (xl(2)-xl(1)) + xl(1);
    yt = 0.75 * (yl(2)-yl(1)) + yl(1);

    set(gca, 'xscale', 'log');
    set(gca, 'yscale', 'log');
    grid on;
    if J==1
        legend('$f_{0,0}$', '$b$', ...
            'FontSize', 12, ...
            'interpreter', 'latex');
    end
    
    % caption = sprintf('$f_{0,0} = %f  x + %f$ \n $b = %f  x + %f$', P_f00(1), P_f00(2), P_bisp(1), P_bisp(2));
    % text(xt, yt, caption, 'FontSize', 12, 'Color', 'black', 'FontWeight', 'bold', ...
    %     'interpreter', 'latex');

    hold off;
end