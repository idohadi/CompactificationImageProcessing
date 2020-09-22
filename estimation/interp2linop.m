function P = interp2linop(td, xsupp, ysupp, im_w, im_h, back_proj)
% TODO: docs in the new convention and update the variable names convention
% and the code base

% Returns a linear interpolation operator P satisfying
%   im_td = P I
% where 
%   I       -   a column vector; the original image, column-major form
%   im_td   -   the interpolated values of the image at the t-design points
%               inside [xsupp(1), xsupp(2)] x [ysupp(1), ysupp(2)] .
% 
% Input arguments
%   td              -   list of spherical t-design points, a point per row, (x,y,z)
%   xsupp, ysupp    -   the image is assumed to be supported in
%                           [xsupp(1), xsupp(2)] x [ysupp(1), ysupp(2)]
%   im_h            -   image height; number of pixels
%   im_w            -   image width; number of pixels
%   back_proj       -   function handle to a projection of points on the 
%                       sphere onto R^2.
%                       Must be expressed in sphe`rical coordinates.
%                       Must support vector input/output:
%                           [r_R2, phi_R2] = back_proj(theta, phi)
%                       Angle convention: 
%                           0<=theta<=pi;
%                           0<=phi<=2pi;
%                           0<=phi_R2<=2pi.
% 
% Output arguments
%   P               -   A matrix satisfying
%                           im_td = P I
%                       as explained above
% 

% Defaults
if nargin<=1
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
x_R2_td = x_R2_td(td_g);
y_R2_td = y_R2_td(td_g);


% Creating a grid
x_grid = linspace(xsupp(1), xsupp(2), im_w); 
y_grid = linspace(ysupp(1), ysupp(2), im_h);
[X, Y] = meshgrid(x_grid, y_grid);

extval = 0; % Value for extrapulation outside x, y

P = zeros(length(x_R2_td), im_w*im_h);
% c = 0;
for w=1:im_w
    S = zeros(length(x_R2_td), im_h);
    parfor h=1:im_h
%         c = (w-1)*im_h + h;
        im_delta = zeros(im_h, im_w);
        im_delta(h, w) = 1;
        S(:, h) = interp2(X, Y, im_delta, x_R2_td, y_R2_td, 'cubic', extval);
    end
    P(:, (w-1)*im_h+1:w*im_h) = S;
    disp(['Completed column ', num2str(w), ' of ', num2str(im_w), '.']);
end

