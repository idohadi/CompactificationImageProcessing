function K = buildK(L, U)
% TODO: docs
% Depends on buildK_mex

UUH = full(U*U');
UUT = full(U*U.');

K = buildK_mex(UUH, UUT, L);
