function [B_r, B_z] = surface_integral_Jungpyo(R, Z, shape, M, varargin)
%{
Calculate the virtual casing principle's surface integral

    1/(4 pi) * \int_{\Gamma} cross( cross(n(y), B(y)), a(x,y) ) d\Gamma(y)

where a(x,y) = (x - y) / || x - y ||^3.
We use the 2D alternating trapezoidal rule in the parameter domain

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

OUTPUTS:
    B_r     (scalar) radial component of the magnetic field at (R,Z)
    B_z     (scalar) vertical """"""""""""""""""""""""""""""""""""""

Evan Toler, 2022
%}

% parse inputs
is_on_surface = false;  if nargin>4, is_on_surface  = varargin{1}; end
t0 = 0;                 if nargin>5, t0             = varargin{2}; end

r = @(t) sqrt( shape.R0^2 + 2*shape.a*shape.R0*cos(t) );
z = @(t) shape.kappa*shape.a*shape.R0 ./ r(t) .* sin(t);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute integrand terms
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dpsi_dr = @(t) shape.kappa*shape.FB/(2*shape.R0^3*shape.q0) .* r(t) ...
    .* ( r(t).^2-shape.R0^2 + 2*z(t).^2./shape.kappa^2 );
dpsi_dz = @(t) shape.a*shape.FB/(shape.R0^2*shape.q0) .* r(t) .* sin(t);

dr_dt = @(t) -shape.a*shape.R0*sin(t) ./ r(t);
dz_dt = @(t) shape.kappa*shape.a*shape.R0 * (r(t).*cos(t) - sin(t).*dr_dt(t)) ./ (r(t)).^2;

% dpsi_dn_bracket = @(t) dpsi_dz(t).*dr_dt(t) - dpsi_dr(t).*dz_dt(t);
% 
alpha = @(t) R^2 + r(t).^2 + (Z-z(t)).^2;
beta  = @(t) 2*R*r(t);
% k2 = @(t) 2*beta(t) ./ ( alpha(t) + beta(t) ); % k squared

Jac = @(t) r(t) .* sqrt(dr_dt(t).^2 + dz_dt(t).^2);
normal = @(t, theta) r(t)./ Jac(t) ...
    .* [dz_dt(t).*cos(theta); 
        dz_dt(t).*sin(theta); 
        -dr_dt(t)];
B_assumed = @(t,theta) 1./r(t) ...
    .* [-dpsi_dz(t).*cos(theta); 
        -dpsi_dz(t).*sin(theta); 
        dpsi_dr(t)];
kernel = @(t,theta) 1./ (alpha(t) - beta(t).*cos(theta)).^(3/2) ...
    .* [R - r(t).*cos(theta); 
        -r(t).*sin(theta);
        (Z - z(t))];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set up integrand and integrate
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% cross product method
h = pi/M;
tt = h*(-M+1:M)';
tt = tt + t0;  % center at singularity
tt = tt - h/2; % straddle the singularity

I = @(t,theta) Jac(t)./(4*pi) * cross(cross(normal(t,theta), B_assumed(t,theta)), kernel(t,theta));
B = zeros(3, 1);
% t = t';
for j = 1:2*M % loop over theta values
    for k = 1:2*M
        B = B + h^2 * I(tt(j), tt(k));
    end
end
B_r = B(1); B_z = B(3);

% % analytic expressions method
% 
% % radial component
% I_r = @(t,theta) 1/(4*pi) .* dpsi_dn_bracket(t) ...
%     .* (Z - z(t)) .* cos(theta) ...
%     ./ (alpha(t) - beta(t).*cos(theta)).^(3/2);
% 
% % z component
% I_z = @(t,theta) 1/(4*pi) .* dpsi_dn_bracket(t) ...
%     .* (r(t) - R.*cos(theta)) ...
%     ./ (alpha(t) - beta(t).*cos(theta)).^(3/2);
% 
% [T, THETA] = meshgrid(tt);
% B_r = h^2 * sum(sum( I_r(T, THETA) ));
% B_z = h^2 * sum(sum( I_z(T, THETA) ));

% add "jump" term when the evaluation point is inside/on surface
if is_on_surface
    dpsi_dr0 = shape.kappa*shape.FB/(2*shape.R0^3*shape.q0) .* R ...
        .* ( R.^2-shape.R0^2 + 2*Z.^2./shape.kappa^2 );
    dpsi_dz0 = shape.a*shape.FB/(shape.R0^2*shape.q0) .* R .* sin(t0);
    B_r = B_r + 0.5 / R * (-dpsi_dz0);
    B_z = B_z + 0.5 / R * dpsi_dr0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot integrand - analytic method only
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% figure;
% surf(T, THETA, I_z(T,THETA))
% xlabel t
% ylabel theta
% zlabel 'B_z Integrand'
% a=.2; xlim([-a,a]); ylim([-a,a])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% check agreement with 1D integral from explicit theta integral
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% omega = (pi - epsi) / 2;
% 
% % I_r_1d = @(t) (1/4*pi) .* dpsi_dn_bracket(t) .* (Z - z(t)) ...
% %     .* (-2)./( alpha(t) + beta(t) ).^(3/2) ...
% %     .* ( ...
% %         2./(1-k2(t)).*ellipticE(omega, k2(t)) ...
% %         - 2.*k2(t)./(1-k2(t)).*sin(omega).*cos(omega)./sqrt(1 - k2(t).*sin(omega).^2) ...
% %         - 4./(k2(t).*(1-k2(t))) .* ellipticE(omega, k2(t)) ...
% %         + 4./k2(t) .* ellipticF(omega, k2(t)) ...
% %         + 4./(1-k2(t)).*sin(omega).*cos(omega)./sqrt(1 - k2(t).*sin(omega).^2) ...
% %         );
% 
% I_r_1d = @(t) (1/4*pi) .* dpsi_dn_bracket(t) .* (Z - z(t)) ...
%     .* (-4)./( alpha(t) + beta(t) ).^(3/2) ...
%     .* ( ...
%         (2-k2(t))./(1-k2(t)).*sin(omega).*cos(omega)./sqrt(1 - k2(t).*sin(omega).^2) ...
%         - (2-k2(t))./(k2(t).*(1-k2(t))) .* ellipticE(omega, k2(t)) ...
%         + 2./k2(t) .* ellipticF(omega, k2(t)) ...
%         );
% 
% B_r_1d = integral(I_r_1d, epsi, 2*pi-epsi ...
%     );
% 
% zterm_1 = @(t) 4.* r(t) ./ (alpha(t) + beta(t)).^(3/2) ...
%     .*  ( ...
%         1./(1-k2(t)) .* ellipticE(omega, k2(t)) ...
%         - k2(t)./(1-k2(t)).*sin(omega).*cos(omega)./sqrt(1 - k2(t).*sin(omega).^2) ...
%         );
% % zterm_2 = @(t) 2*R./ (alpha(t) + beta(t)).^(3/2) ...
% %     .*  ( ...
% %         2./(1-k2(t)).*ellipticE(omega, k2(t)) ...
% %         - 2.*k2(t)./(1-k2(t)).*sin(omega).*cos(omega)./sqrt(1 - k2(t).*sin(omega).^2) ...
% %         - 4./(k2(t).*(1-k2(t))) .* ellipticE(omega, k2(t)) ...
% %         + 4./k2(t) .* ellipticF(omega, k2(t)) ...
% %         + 4./(1-k2(t)).*sin(omega).*cos(omega)./sqrt(1 - k2(t).*sin(omega).^2) ...
% %         ); % copy/paste from the cos/(...)^(3/2) integral
% zterm_2 = @(t) 4*R./ (alpha(t) + beta(t)).^(3/2) ...
%     .*  ( ...
%         (2-k2(t))./(1-k2(t)).*sin(omega).*cos(omega)./sqrt(1 - k2(t).*sin(omega).^2) ...
%         - (2-k2(t))./k2(t)./(1-k2(t)).*ellipticE(omega, k2(t)) ...
%         + 2./k2(t) .* ellipticF(omega, k2(t)) ...
%         );
%    
% I_z_1d = @(t) (1/4*pi) .* dpsi_dn_bracket(t) ...
%     .* (zterm_1(t) + zterm_2(t));
% B_z_1d = integral(I_z_1d, epsi, 2*pi-epsi ...
%     );
% 
% % add jump
% % TO DO AFTER UN-COMMENTING
%
% fprintf('|B_r - B_r_1d| = %1.12e\n', abs(B_r-B_r_1d))
% fprintf('|B_z - B_z_1d| = %1.12e\n', abs(B_z-B_z_1d))
% fprintf('Zakharov line integral:\n')
% fprintf('[B_r=%1.15e, B_z=%1.15e] at (r=%1.1f, z=%1.1f)\n', B_r, B_z, R, Z)

% % plot 1D integrand
% figure(1); clf;
% tt = linspace(epsi, 2*pi-epsi, 1e2);
% plot(tt, I_r_1d(tt)) % plot integrand
% xlabel $t$; ylabel '$B_r$ Integrand'; title 'Line Integral on surface'
% xlim([0, 2*pi])
% 
% figure(2); clf;
% plot(tt, I_z_1d(tt)) % plot integrand
% xlabel $t$; ylabel '$B_z$ Integrand'; title 'Line Integral on surface'
% xlim([0, 2*pi])
% 
% % plot singularity near t=0
% figure(3); clf;
% tt = logspace(-1,-10,10);
% loglog(tt, abs(I_r_1d(tt)), 'x-', tt, abs(I_z_1d(tt)), 'o-') % plot integrand near 0
% axis tight
% hold on
% loglog(tt, abs(log(tt)), 'k--')
% xlabel $t$
% ylabel 'Integrands'
% legend('B_r integrand', 'B_z integrand', '-log(t)')
return
