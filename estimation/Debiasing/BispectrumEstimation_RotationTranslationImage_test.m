%% Parameter initialization
% Band-limit
L = 5;

% Load t-design
t = 70;
td = loadtd('sf070.02522');

% Support (boundary) of image
xsupp = [-0.5, 0.5];
ysupp = xsupp;

% Image size; in pixels
im_w = 101;     % width
im_h = im_w;    % height

% Back-projectino from the sphere into R^2; in spherical coordinates
a = 2;
back_proj = @(theta, phi) deal(theta/a, phi);

% U matrix, containing debiasing coefficients
U = debiasmatrix(td, L, xsupp, ysupp, im_w, im_h, back_proj);

% Loading Clebsch-Gorda coefficients
CGs = build_CGs_vec(L);

%% Generate a random image
im = rand(im_h, im_w);

% add padding
w_pad = round(im_w/2);
h_pad = round(im_h/2);

im = padarray(im,[h_pad, w_pad],0,'both');


%% Slope test parameters

sigma = [1];
sigma = sigma(:);

N = 100:5:200;
N = N(:);

trial_num = 10;

%% The slope test
im_shc = image2shc(im, L, td, back_proj, xsupp, ysupp);
b_im = bisp_vec(im_shc, L, CGs);

results = zeros(length(N), length(sigma));
t1_all = tic();
for s=1:length(sigma)
    
    disp(strcat(['Sigma = ', num2str(sigma(s))]));
    t1_sigma = tic();
    for n=1:length(N)
        
        disp(strcat(['Exeriment number ', num2str(n), ' of ', num2str(length(N)), '. N = ', num2str(N(n))]));
        
        t1_N = tic();
        
        trial_res = [];
        for k=1:trial_num
            
            disp(strcat(['Trial number ', num2str(k), ' of ', num2str(trial_num)]));
            
            t1 = tic();
            
            im_noisy = noisyRotatedTranslatedImage(im, N(n), sigma(s));
            im_shc = image2shc(im_noisy, L, td, back_proj, xsupp, ysupp);
            b = estimbisp(im_shc, L, CGs, U, sigma(s));
            trial_res(k) = max(abs(b - b_im))/max(abs(b_im));
            
            t2 = toc(t1);
            disp(strcat(['Trial duration: ', num2str(t2)]));
        end
        
        results(n, s) = mean(trial_res);
        
        t2_N = toc(t1_N);
        disp(strcat(['Experiment duration: ', num2str(t2_N)]));
    end
    
    t2_sigma = toc(t1_sigma);
    disp(strcat(['Overall for sigma: ', num2str(t2_sigma)]));
end
t2_all = toc(t1_all);
disp(strcat(['Overall time: ', num2str(t2_all)]));

%% Plotting the slope test
figure;
for J=1:1
    
    subplot(1, 1, J);
    
    scatter(N, results(:, J), 12, 'filled', 'blue'); 

    P_bisp = polyfit(log(N), log(results(:, J)), 1);

    xlabel('$N$, sample size', ...
        'FontSize', 14, ...
        'interpreter', 'latex');
    ylabel({'max-norm relative', 'to ground truth'}, ...
        'FontSize', 14, ...
        'interpreter', 'latex');
    caption = strcat(['$\sigma^2 = ', num2str(sigma(J)), '$']);
    title({caption, ['$b$ slope: ', num2str(P_bisp(1))]}, ...
        'Fontsize', 13, ...
        'interpreter', 'latex');

    set(gca, 'xscale', 'log');
    set(gca, 'yscale', 'log');
    grid on;
end
