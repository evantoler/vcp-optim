function Tf = periodic_KR(h, ff, order)
%{
Compute \int_{-b}^b f(t) dt using a trapezoidal
rule with Kapur-Rokhlin corrections when f is 2b-periodic
and has a logarithmic singularity at the origin

INPUTS:
    h       (scalar) constant grid spacing (=2b/2M)
    ff      (2M) integrand evaluated at equispaced values; that is,
                ff(j) = f( -b+j*h ).
                The singularity should occur at index M
    order   (in {2,6,10}) The KR convergence order parameter.

OUTPUT:
    Tf      (scalar) estimate of the integral in the preamble

Evan Toler, 2021.
%}

% singularity corrections near the origin; don't need endpoint corrections
gamma = KR_weights(order);

M = length(ff)/2;

Tf = h * sum( ff([1:M-1, M+1:end]) ); % standard "punctured" trapezoidal

% singularity corrections; points within "order" of the singularity
% are corrected twice by symmetry
if order>0
    idx = [-order:-1, 1:order];
    Tf = Tf + h * dot(gamma + flip(gamma), ff( M+idx ));
end
