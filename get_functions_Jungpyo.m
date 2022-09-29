function [rr, zz, dpsi_dr, dpsi_dz, dr_dt, dz_dt, dpsi_dn_bracket, ...
    alpha, beta, k2, d2r_dt2, d2z_dt2] = get_functions_Jungpyo(shape, t0, tt)
%{
Generate points and derivatives wrt the generating curve parameterization.
Also compute point-to-point quantities associated with another point on the
generating curve.

INPUTS:
    shape   (struct) a shape object representing a geometry
                parameterization from Lee & Cerfon.
    t0      (scalar) parameter of the associated surface point for
                point-to-point quantities
    tt      (N) vector of parameters at which to evaluate the points and
                derivatives; point-to-point quantities are computed between
                t0 and tt(j) for each j in 1:N

OUTPUTS:
    rr, zz              (N) points on the generating curve; jth element is
                            at tt(j).
    dpsi_dr, dpsi_dz    (N) derivative of flux function wrt r, z evaluated
                            at rr(1:N), zz(1:N)
    dr_dt, dz_dt        (N) derivative of points wrt paramter t at tt(1:N)
    dpsi_dn_bracket     (N) normal derivative of flux function evaluated at
                            tt(1:N)
    alpha, beta, k2     (N) point-to-point quantities b/t tt(1:N) and t0;
                            example: beta(j) = 2*r(t0)*r(tt(j))
    d2r_dt2, d2z_dt2    (N) analogs of dr_dt, dz_dt for 2nd derivatives

Evan Toler, 2022
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% discretization of points on the generating curve
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rr = sqrt( shape.R0^2 + 2*shape.a*shape.R0*cos(tt) );
zz = shape.kappa*shape.a*shape.R0 ./ rr .* sin(tt);

if nargout < 3, return; end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1st derivatives
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% observation point
R = sqrt( shape.R0^2 + 2*shape.a*shape.R0*cos(t0) );
Z = shape.kappa*shape.a*shape.R0 ./ R .* sin(t0);

dpsi_dr = shape.kappa*shape.FB/(2*shape.R0^3*shape.q0) .* rr ...
    .* ( rr.^2 - shape.R0^2 + 2*zz.^2./shape.kappa^2 );
dpsi_dz = shape.FB/(shape.kappa*shape.R0^3*shape.q0) .* rr.^2 .* zz;
% alternative formula by plugging in zz directly
% dpsi_dz = shape.a*shape.FB/(shape.R0^2*shape.q0) .* rr .* sin(tt);

dr_dt = -shape.a*shape.R0*sin(tt) ./ rr;
dz_dt = shape.kappa*shape.a*shape.R0 * (cos(tt)./rr - sin(tt)./rr.^2.*dr_dt);

if nargout < 7, return; end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1D Integration / Elliptic Integral Numbers
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dpsi_dn_bracket = dpsi_dz.*dr_dt - dpsi_dr.*dz_dt;

alpha = R^2 + rr.^2 + (Z-zz).^2;
beta = 2*R.*rr;
k2 = 2.*beta ./ ( alpha + beta ); % k squared

if nargout < 11, return; end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2nd derivatives (to check analysis / debugging)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

d2r_dt2 = -shape.a*shape.R0 * (cos(tt) ./ rr - sin(tt)./rr.^2.*dr_dt);
d2z_dt2 = shape.kappa*shape.a*shape.R0 * ( -sin(tt)./rr - cos(tt)./rr.^2.*dr_dt ...
    - cos(tt)./rr.^2.*dr_dt + 2*sin(tt)./rr.^3.*dr_dt.^2 - sin(tt)./rr.^2.*d2r_dt2 );