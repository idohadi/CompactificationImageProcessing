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

%% Obtain t-design points within the support
% Back-project t-design onto R^2
[theta_sp_td, phi_sp_td] = sphcart2sph(td(:, 1), td(:, 2), td(:, 3));
[r_R2_td, phi_R2_td] = back_proj(theta_sp_td, phi_sp_td);
[x_R2_td, y_R2_td] = pol2cartm(r_R2_td, phi_R2_td);

% Isolating t-design points falling within the support, defined by xsupp and ysupp
td_g = x_R2_td>=xsupp(1) & x_R2_td<=xsupp(2) & y_R2_td>=ysupp(1) & y_R2_td<=ysupp(2);
x_R2_td = x_R2_td(td_g);
y_R2_td = y_R2_td(td_g);

%% Calculate P matrix
P = interp2linop(td, xsupp, ysupp, im_w, im_h, back_proj);
P = sparse(P);

%% Tests for interp2linop

% random image
im = rand(im_h, im_w);
% im = zeros(im_h, im_w);
% im(randi(im_h), randi(im_w)) = 1;
extval = 0; % Value for extrapulation outside support

im_vals_interp = interp2(X, Y, im, x_R2_td, y_R2_td, 'cubic', extval);

im_vals_P = P * im(:);

max(abs(im_vals_interp - im_vals_P))
% Difference is O(10^(-16))


%% Calculate the U matrix
U = debiasmatrix(td, L, xsupp, ysupp, im_w, im_h, back_proj);

%% Tests for debiasmatrix

im = rand(im_h, im_w);

im_shc = image2shc(im, L, td, back_proj, xsupp, ysupp);

im_shc_U = U * im(:);

max(abs(im_shc - im_shc_U))
% Difference is  O(10^(-16))
