function [B_r, B_z] = Zakharov_line_integral_Jungpyo(R, Z, shape, M, varargin)
%{
Calculate the virtual casing principle's surface integral

    1/(4 pi) * \int_{\Gamma} cross( cross(n(y), B(y)), a(x,y) ) d\Gamma(y)

where a(x,y) = (x - y) / || x - y ||^3.

INPUTS:
    R,Z             (scalars) target point location
    shape           (struct) geometry parameters with fields as outlined 
                        in Lee & Cerfon:
                        .kappa
                        .a
                        .R0
                        .FB
                        .q0
    M               (int >0) integrate using 2M quadrature points
    is_on_surface   (boolean) true if (R,Z) is on the surface; else false
    t0              (scalar) if (is_on_surface), this is the value of
                        the parameter such that (R,Z) = (r(t0), z(t0))
                        default 0
    alternating     (boolean) true if (is_on_surface) and the user 
                        wishes to use alternating trapeoidal rule
                        quadrature
                        default false
    order           (in {2,6,10}) Kapur-Rokhlin parameter; see KR_weights.m
                        default 10

OUTPUTS:
    B_r     (scalar) radial component of the magnetic field at (R,Z)
    B_z     (scalar) vertical """"""""""""""""""""""""""""""""""""""

Evan Toler, 2022
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% parameterize the surface
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

h = pi/M;
tt = h.*(-M+1 : M)'; %(-pi+h), ..., pi

% default values
is_on_surface = false;  if nargin>4, is_on_surface  = varargin{1}; end
t0 = 0;                 if nargin>5, t0             = varargin{2}; end
alternating = false;    if nargin>6, alternating    = varargin{3}; end
order = 10;             if nargin>7, order          = varargin{4}; end

% if user doesn't know if the point is on the surface...
if nargin <= 4
    [rr,zz] = get_functions_Jungpyo(shape, t0, tt); % default t0
    dists = (R-rr).^2 + (Z-zz).^2;
    is_on_surface = any(dists < 1e-3);
    [~, idx] = min(dists);
    t0 = tt(idx);
end

% parse inputs
if is_on_surface

    % shift parameterization so singularity (t0) is at middle index
    tt = tt + t0;

    if alternating
        % fprintf('Line integral: on surface. Using an alternating (staggered/midpoint) quad rule.\n')
        tt = tt - h/2; % stagger
    else
        % fprintf('Line integral: on surface. Using a Kapur-Rokhlin quad rule.\n')
    end

end

% integration points
[rr,zz] = get_functions_Jungpyo(shape, t0, tt);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute integrand terms
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dpsi_dr = shape.kappa*shape.FB/(2*shape.R0^3*shape.q0) .* rr ...
    .* ( rr.^2-shape.R0^2 + 2*zz.^2./shape.kappa^2 );
dpsi_dz = shape.a*shape.FB/(shape.R0^2*shape.q0) .* rr .* sin(tt);

dr_dt = -shape.a*shape.R0*sin(tt) ./ rr;
dz_dt = shape.kappa*shape.a*shape.R0 * (cos(tt)./rr - sin(tt)./rr.^2.*dr_dt);

dpsi_dn_bracket = dpsi_dz.*dr_dt - dpsi_dr.*dz_dt;

k2 = 4*R.*rr ./ ( (R+rr).^2 + (Z-zz).^2 ); % k squared

% radial component of integrand
Ecoef = R^2 + rr.^2 + (Z-zz).^2;
Ecoef = Ecoef ./ ( (R-rr).^2 + (Z-zz).^2 );
B_r = 1/(2*pi) ./ rr .* dpsi_dn_bracket .* (Z-zz)./R;
B_r = B_r ./  sqrt( (R+rr).^2 + (Z-zz).^2 );
B_r = B_r .* ( -ellipticK(k2) + Ecoef.*ellipticE(k2) );

% z component of integrand
Ecoef = rr.^2 - R^2 - (Z-zz).^2;
Ecoef = Ecoef ./ ( (rr-R).^2 + (Z-zz).^2 );
B_z = 1/(2*pi) ./ rr .* dpsi_dn_bracket;
B_z = B_z ./ sqrt( (R+rr).^2 + (Z-zz).^2 );
B_z = B_z .* ( ellipticK(k2) + Ecoef.*ellipticE(k2) );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot integrand
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% if is_on_surface
%     figure;
%     clf;
%     plot(tt([1:M-1. M+1:end]), B_z([1:M-1. M+1:end])) % plot integrand
%     xlabel $t$; ylabel $f(t)$;
%     title(sprintf('Integrand, $B_z = 2 \\cdot \\frac 1{4\\pi} \\int_{-\\pi}^\\pi f(t) \\, dt$;   $(R, Z)=(%1.1e, %1.1e)$', R,Z))
%     xlim([-pi, pi])
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute line integral
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if is_on_surface % corrected trapezoidal rule (alternating or KR)
    
    if alternating % alternating trapezoidal
        B_r = h*sum(B_r);
        B_z = h*sum(B_z);
    else % KR
        B_r = periodic_KR(h, B_r, order);
        B_z = periodic_KR(h, B_z, order);
    end
    
    % add the jump
    dpsi_dz = shape.a*shape.FB/(shape.R0^2*shape.q0) .* R .* sin(t0); % evaluated at t = t0
    dpsi_dr = shape.kappa*shape.FB/(2*shape.R0^3*shape.q0) * R ...
        .* ( R^2-shape.R0^2 + 2*Z^2/shape.kappa^2 );
    B_r = B_r + 0.5 / R * (-dpsi_dz);
    B_z = B_z + 0.5 / R * dpsi_dr;
    
else % periodic trapezoidal rule
    
    B_r = h*sum(B_r);
    B_z = h*sum(B_z);

end
