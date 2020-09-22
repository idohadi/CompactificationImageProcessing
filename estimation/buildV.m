function V = buildV(L, U)
% TODO: docs

V = buildV_mex(U*U', L);
V = real(V);
