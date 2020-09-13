%% Parameter initialization
% Band-limit
L = 8;

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

% Creating a grid
x_grid = linspace(xsupp(1), xsupp(2), im_w); 
y_grid = linspace(ysupp(1), ysupp(2), im_h);
[X, Y] = meshgrid(x_grid, y_grid);

% U matrix, containing debiasing coefficients
U = debiasmatrix(td, L, xsupp, ysupp, im_w, im_h, back_proj);

% Loading Clebsch-Gorda coefficients
CGs = build_CGs_vec(L);

%% Dry test, to see that the funciton works
% The image has no noise, sigma = 0

im_shc = rand((L+1)^2, 1);

sigma = 0;
b_est = estimbisp(im_shc, L, CGs, U, sigma);
b_real = bisp_vec(im_shc, L, CGs);


max(abs(b_est-b_real)/max(abs(b_real)))
