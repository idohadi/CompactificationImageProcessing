function U = buildU(L, td, xsupp, ysupp, im_w, im_h, back_proj)
% TODO: docs and update the variable names to my current convention and
% also update the notations

% Back-project t-design onto R^2
[theta_sp_td, phi_sp_td] = sphcart2sph(td(:, 1), td(:, 2), td(:, 3));
[r_R2_td, phi_R2_td] = back_proj(theta_sp_td, phi_sp_td);
[x_R2_td, y_R2_td] = pol2cartm(r_R2_td, phi_R2_td);

% Isolating t-design points falling within the support, defined by xsupp and ysupp
td_g = x_R2_td>=xsupp(1) & x_R2_td<=xsupp(2) & y_R2_td>=ysupp(1) & y_R2_td<=ysupp(2);
theta_sp_td = theta_sp_td(td_g);
phi_sp_td = phi_sp_td(td_g);

Yt = evalYt(theta_sp_td, phi_sp_td, L);

U = 4*pi*conj(Yt)*interp2linop(td, xsupp, ysupp, im_w, im_h, back_proj)/size(td, 1);

U = sparse(U);
