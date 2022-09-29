%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Verify derivatives of the ring charge potential.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear

syms r r0 z z0

% fundamental quantities
alpha = r^2 + r0^2 + (z-z0)^2;
beta = 2 * r * r0;
k2 = 2*beta / (alpha+beta);

% ring charge potential function
u = 4 * r0 / sqrt(alpha + beta) * ellipticK(k2);

% verify that u is harmonic
% (It should be, as a linear superposition of harmonic functions.)
Lapu = r * diff(u,r);
Lapu = diff(Lapu,r)/r;
Lapu = Lapu + diff(u,z,z);

Lapu = simplify(Lapu);
pretty(Lapu)

% print expressions for gradient components
% dudr = diff(u,r); dudr = simplify(dudr);
% dudz = diff(u,z); dudz = simplify(dudz);
% pretty(dudr)