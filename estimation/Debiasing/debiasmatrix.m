function U = debiasmatrix(td, L, xsupp, ysupp, im_w, im_h, back_proj)
% Constructs the debiasing coefficients matrix U.
% 
% Input arguments
%   td              -   list of spherical t-design points, a point per row, (x,y,z)
%   L               -   band-limit of SHCs that can be calculated using U
%   xsupp, ysupp    -   the image is assumed to be supported in
%                           [xsupp(1), xsupp(2)] x [ysupp(1), ysupp(2)]
%   im_h            -   image height; number of pixels
%   im_w            -   image width; number of pixels
%   
% 
% Output arguments
%   U       -   the debiasing matrix
% 
% NOTE:
%   This assumes the noise is in the image space.
% 

% Defaults
if nargin<=2
    xsupp = [-0.5, 0.5];
    ysupp = [-0.5, 0.5];
    
    im_w = 101;
    im_h = 101;
    
    a = 2;
    back_proj = @(theta, phi) deal(theta/a, phi);
end

% Back-project t-design onto R^2
[theta_sp_td, phi_sp_td] = sphcart2sph(td(:, 1), td(:, 2), td(:, 3));
[r_R2_td, phi_R2_td] = back_proj(theta_sp_td, phi_sp_td);
[x_R2_td, y_R2_td] = pol2cartm(r_R2_td, phi_R2_td);

% Isolating t-design points falling within the support, defined by xsupp and ysupp
td_g = x_R2_td>=xsupp(1) & x_R2_td<=xsupp(2) & y_R2_td>=ysupp(1) & y_R2_td<=ysupp(2);
theta_sp_td = theta_sp_td(td_g);
phi_sp_td = phi_sp_td(td_g);

% Calculate Y_(t') matrix
Yttag = evalYt(theta_sp_td, phi_sp_td, L);

% Calculates the linear operator of the inerpolation
P = interp2linop(td, xsupp, ysupp, im_w, im_h, back_proj);

U = 4*pi/length(td) * conj(Yttag) * P;
