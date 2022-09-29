%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% examine the integrands involved with a principal value integral
% using the Zakharov line integral method
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 0. Set up shape and observation point
clear; close all; clc

M = 100;            % number of discretization points = 2*M
epsi = 1e-2;        % p.v. parameter
version = 'normal'; % 'log' to see behavior near boundary

% set up shape
shape.kappa = 1;
shape.a = 1/3;
shape.R0 = 1;
shape.FB = 1;
shape.q0 = 1;
Rmin = sqrt( shape.R0^2 - 2*shape.a*shape.R0 );
Rmax = sqrt( shape.R0^2 + 2*shape.a*shape.R0 );

% set observation point -- at t=t0
t0 = 1;
R = sqrt( shape.R0^2 + 2*shape.a*shape.R0*cos(t0) );
Z = shape.kappa*shape.a*shape.R0 ./ R .* sin(t0);

% principal value parameter -- excising a square in the parameter domain
omega = (pi-epsi)/2;

%% 1. Parameterize the Surface and get Important Quantities

switch version
    case 'log'
        tt = logspace(-10, 0, 2*M)' + t0;
    otherwise
        h = pi/M;
        tt = h.*(-M+1 : M)' + t0; %(-pi+h), ..., pi shifted by t0
end

[rr, ~, ~, ~, ~, dz_dt, dpsi_dn_bracket, ...
    alpha, beta, k2, d2r_dt2, d2z_dt2] = get_functions_Jungpyo(shape, t0, tt);

%% 2. Define integrands in respective regions

C = 1/(4*pi);

% Determine where parameters are near the singularity
idx_small = ( abs(tt - t0) <= epsi ); 

% unadjusted
% Iz = C * dpsi_dn_bracket * 4./(alpha+beta).^(3/2);

% for t away from singularity at t0
fz_far = C * dpsi_dn_bracket * 4./(alpha+beta).^(3/2) ...
    .* ( ...
    ( (rr - R)./(1-k2) - 2*R./k2 ) .* ellipticE(k2) ...
    + 2*R./k2 .* ellipticK(k2) ...
    );
fz_far( idx_small) = 0;

% for t near the singularity at t0
fz_near = C * dpsi_dn_bracket * 4./(alpha+beta).^(3/2) ...
    .* ( ...
    ( (rr - R)./(1-k2) - 2*R./k2 ) .* ellipticE(omega, k2) ...
    + 2*R./k2 .* ellipticF(omega, k2) ...
    + ( (R-rr).*k2./(1-k2) + 2*R ) ...
        .* sin(omega)*cos(omega)./sqrt(1-k2.*(sin(omega)^2)) ...
    );
fz_near(~idx_small) = 0;
fz_near(M) = C * dpsi_dn_bracket(M) * 4./(alpha(M)+beta(M)).^(3/2) ...
    .* ( ...
    ( 2*R^2*d2r_dt2(M)/dz_dt(M)^2 - 2*R./k2(M) ) .* ellipticE(omega, k2(M)) ...
    + 2*R./k2(M) .* ellipticF(omega, k2(M)) ...
    + ( -2*R^2*d2r_dt2(M)/dz_dt(M)^2 .* k2(M) + 2*R ) ...
        .* sin(omega)*cos(omega)./sqrt(1-k2(M).*(sin(omega)^2)) ...
    ); % limiting value as t -> t0; removable singularity

%% 3. Plot Integrands

figure(1); clf
switch version
    case 'log'
        legendstr = {'$|t-t_0|>\varepsilon$', '$|t-t_0|\le \varepsilon$', '$1/|t-t_0|$'};
        loglog(tt-t0, abs(fz_far), tt-t0, abs(fz_near), tt-t0, abs(1./(tt-t0)),'k--')
        xlabel '$t - t_0$'
    otherwise
        legendstr = {'$|t-t_0|>\varepsilon$', '$|t-t_0|\le \varepsilon$'};
        plot(tt, fz_far, tt, fz_near)
        xlabel '$t$'
        xlim([-pi, pi] + t0)
end

ylabel('$f_z(t; \epsilon)$','interpreter','latex')
title(sprintf('$\\varepsilon = %1.1e$', epsi),'interpreter','latex')
legend(legendstr, 'interpreter', 'latex')

%% 4a. Built-in integration -- making function handles

clear tt
rr = @(t) sqrt( shape.R0^2 + 2*shape.a*shape.R0*cos(t) );
zz = @(t) shape.kappa*shape.a*shape.R0 ./ rr(t) .* sin(t);

dpsi_dr = @(t) shape.kappa*shape.FB/(2*shape.R0^3*shape.q0) .* rr(t) ...
    .* ( rr(t).^2-shape.R0^2 + 2*zz(t).^2./shape.kappa^2 );
dpsi_dz = @(t) shape.a*shape.FB/(shape.R0^2*shape.q0) .* rr(t) .* sin(t);

dr_dt = @(t) -shape.a*shape.R0*sin(t) ./ rr(t);
dz_dt = @(t) shape.kappa*shape.a*shape.R0 * (cos(t)./rr(t) - sin(t)./rr(t).^2.*dr_dt(t));

dpsi_dn_bracket = @(t) dpsi_dz(t).*dr_dt(t) - dpsi_dr(t).*dz_dt(t);

alpha = @(t) R^2 + rr(t).^2 + (Z-zz(t)).^2;
beta = @(t) 2*R.*rr(t);
k2 = @(t) 2*beta(t) ./ ( alpha(t) + beta(t) ); % k squared

% z-component of integrand
f = @(t,phi) C * dpsi_dn_bracket(t) .* (alpha(t)-beta(t).*cos(phi)).^(-3/2) .* (rr(t) - R*cos(phi));

%% 4b. Built-in integration -- integrate and plot integrand

% pv parameters
nepsi = 20;
epsiepsi = logspace(-10, 0, nepsi)';

% preallocate
[I1I1, I2I2] = deal(zeros(size(epsiepsi)));

for i = nepsi:-1:1
    
    fprintf('iter %d\n',i)
    epsi = epsiepsi(i);

    % built-in integral on most of rectangle
    I1 = integral2(f, t0+epsi, t0+2*pi-epsi, -pi, pi, 'abstol',1e-15, 'reltol',1e-15);
    
    % built-in integral on the thin strip minus rectangle around singularity
    I2 = integral2(f, t0-epsi, t0+epsi, epsi, 2*pi-epsi, 'abstol',1e-15, 'reltol',1e-15);

    I1I1(i) = I1;
    I2I2(i) = I2;
    
    % visualize integrand near the singularity
    [T,PHI] = meshgrid(logspace(-10,0,100)+t0, logspace(-10, 0, 100));
    figure(2); clf
    fz = f(T,PHI);
    surf(log(T-t0), log(PHI), log(abs(fz)))
    xlabel('$\log(t - t_0)$','interpreter','latex')
    ylabel('$\log(\phi)$','interpreter','latex')
    zlabel('$\log( |e_z \cdot \vec F(\theta, t) |)$','interpreter','latex')
    title(sprintf('$\\varepsilon = %1.1e$', epsi),'interpreter', 'latex')
    hold on
    surf(log(T-t0), log(PHI), log(abs( 1./((T-t0).*(PHI-0)) )) )
    shading interp
    legend('2D surface integrand', '$1/[(t-t_0)(\phi-0)]$','interpreter','latex')

    % debugging
    % bar = 2*epsi*max(max(abs(fz)));
    % if i==22, keyboard; end
    % disp( bar ) % upper bound on integral; should -> 0
    
    drawnow;
end

figure(3); clf
loglog(epsiepsi, abs(I1I1), 'x-'); hold on
loglog(epsiepsi, abs(I2I2), 'o-')
loglog(epsiepsi, abs(I1I1 + I2I2), 's-')
legend('Integral over most of $(t,\theta)$ [in abs()]', ...
    'Integral on thin strip $t \in [-\epsilon,\epsilon]$ [in abs()]', ...
    'Sum; total PV integral [in abs()]', ...
    'interpreter', 'latex')
xlabel('$\epsilon$','interpreter','latex')

%% 5. Test alternating trapezoidal rule

epsi = 1e-4; % PV parameter

% "true" integral values
I1 = integral2(f, t0+epsi, t0+2*pi-epsi, -pi, pi, 'abstol',1e-15, 'reltol',1e-15);
I2 = integral2(f, t0-epsi, t0+epsi, epsi, 2*pi-epsi, 'abstol',1e-15, 'reltol',1e-15);

% copy/paste from "define integrands in respective regions"
nM = 40;
MM = round(exp(linspace(log(5), log(500), nM)))';
[I1I1alt,I2I2alt] = deal( zeros(size(MM)) );

for j = 1:nM
    
    M = MM(j);
    h = pi/M;
    tt = h.*(-M+1 : M)' + t0; %(-pi+h), ..., pi shifted by t0
    tt = tt - h/2;

    [rr, zz, dpsi_dr, dpsi_dz, dr_dt, dz_dt, dpsi_dn_bracket, ...
        alpha, beta, k2, d2r_dt2, d2z_dt2] = get_functions_Jungpyo(shape, t0, tt);

    idx_small = ( abs(tt - t0) <= epsi );
    C = 1/(4*pi);
    
    % unadjusted
    % Iz = C * dpsi_dn_bracket * 4./(alpha+beta).^(3/2);
    
    % for t away from singularity at t0
    fz_far = C * dpsi_dn_bracket * 4./(alpha+beta).^(3/2) ...
        .* ( ...
        ( (rr - R)./(1-k2) - 2*R./k2 ) .* ellipticE(k2) ...
        + 2*R./k2 .* ellipticK(k2) ...
        );
    fz_far( idx_small) = 0;
    
    % for t near the singularity at t0
    fz_near = C * dpsi_dn_bracket * 4./(alpha+beta).^(3/2) ...
        .* ( ...
        ( (rr - R)./(1-k2) - 2*R./k2 ) .* ellipticE(omega, k2) ...
        + 2*R./k2 .* ellipticF(omega, k2) ...
        + ( (R-rr).*k2./(1-k2) + 2*R ) ...
        .* sin(omega)*cos(omega)./sqrt(1-k2.*(sin(omega)^2)) ...
        );
    fz_near(~idx_small) = 0;

    % debugging: check limiting value at the singularity
%     fz2(M) = C * dpsi_dn_bracket(M) * 4./(alpha(M)+beta(M)).^(3/2) ...
%         .* ( ...
%         ( 2*R^2*d2r_dt2(M)/dz_dt(M)^2 - 2*R./k2(M) ) .* ellipticE(omega, k2(M)) ...
%         + 2*R./k2(M) .* ellipticF(omega, k2(M)) ...
%         + ( -2*R^2*d2r_dt2(M)/dz_dt(M)^2 .* k2(M) + 2*R ) ...
%         .* sin(omega)*cos(omega)./sqrt(1-k2(M).*(sin(omega)^2)) ...
%         ); % limiting value as t -> 0; removable singularity

    %%%%%%%%%%% integrate %%%%%%%%%%%
    I1I1alt(j) = h * sum(fz_far);
    I2I2alt(j) = h * sum(fz_near);
    
end

figure(4); clf
loglog(MM, abs(I1I1alt - I1) ./ abs(I1), 'x-'); hold on
loglog(MM, abs(I2I2alt - I2) ./ abs(I2), 'o-')
loglog(MM, 1./MM, 'k--')
legend('Error in bulk integral', 'Error in thin strip integral', 'first order')
xlabel('$M$ ($2M$ quadrature points)','interpreter','latex')
ylabel('Alternating trapezoidal rule error')
title(sprintf('$\\varepsilon = %1.1e$', epsi),'interpreter', 'latex')
axis tight

%% 6. Kapur-Rokhlin
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% test Kapur-Rokhlin trapezoidal correction for PV integral
% (not an appropriate quadrature rule -- just for comparison)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% KR convergence order parameter
order = 10;

% PV parameters for approximate checking
epsi = 1e-4;

% "true" integral values
I1 = integral2(f, t0+epsi, t0+2*pi-epsi, -pi, pi, 'abstol',1e-15, 'reltol',1e-15);
I2 = 0;

% convergence scaling in M (2M = # quadrature nodes)
nM = 40;
MM = round(exp(linspace(log(order+1), log(500), nM)))';

% preallocate
[I1I1KR,I2I2KR] = deal( zeros(size(MM)) );

for j = 1:nM
    
    fprintf('Iter %d/%d\n', j, nM)
    
    M = MM(j);
    h = pi/M;
    tt = h.*(-M+1 : M)' + t0; %(-pi+h), ..., pi shifted by t0

    % copy/paste from "define integrands in respective regions"
    
    [rr, zz, dpsi_dr, dpsi_dz, dr_dt, dz_dt, dpsi_dn_bracket, ...
        alpha, beta, k2, d2r_dt2, d2z_dt2] = get_functions_Jungpyo(shape, t0, tt);
    
    idx_small = ( abs(tt - t0) <= epsi );
    C = 1/(4*pi);
    
    % for t away from singularity at t0
    fz_far = C * dpsi_dn_bracket * 4./(alpha+beta).^(3/2) ...
        .* ( ...
        ( (rr - R)./(1-k2) - 2*R./k2 ) .* ellipticE(k2) ...
        + 2*R./k2 .* ellipticK(k2) ...
        );
    fz_far( idx_small) = 0;
    
    % for t near the singularity at t0
    fz_near = C * dpsi_dn_bracket * 4./(alpha+beta).^(3/2) ...
        .* ( ...
        ( (rr - R)./(1-k2) - 2*R./k2 ) .* ellipticE(omega, k2) ...
        + 2*R./k2 .* ellipticF(omega, k2) ...
        + ( (R-rr).*k2./(1-k2) + 2*R ) ...
        .* sin(omega)*cos(omega)./sqrt(1-k2.*(sin(omega)^2)) ...
        );
    fz_near(~idx_small) = 0;
    fz_near(M) = C * dpsi_dn_bracket(M) * 4./(alpha(M)+beta(M)).^(3/2) ...
        .* ( ...
        ( 2*R^2*d2r_dt2(M)/dz_dt(M)^2 - 2*R./k2(M) ) .* ellipticE(omega, k2(M)) ...
        + 2*R./k2(M) .* ellipticF(omega, k2(M)) ...
        + ( -2*R^2*d2r_dt2(M)/dz_dt(M)^2 .* k2(M) + 2*R ) ...
        .* sin(omega)*cos(omega)./sqrt(1-k2(M).*(sin(omega)^2)) ...
        ); % limiting value as t -> t0; removable singularity

    %%%%%%%%%%% integrate %%%%%%%%%%%
    I1I1KR(j) = periodic_KR(h, fz_far, order);
    I2I2KR(j) = periodic_KR(h, fz_near, order);
    
end

figure(5); clf
loglog(MM, abs(I1I1KR - I1) ./ abs(I1), 'x-'); hold on
loglog(MM, abs(I2I2KR - I2), 'o-')
loglog(MM, 1./MM, 'k--')
legend('Error in bulk integral', 'Error in thin strip integral', 'first order')
xlabel('$M$ ($2M$ quadrature points)','interpreter','latex')
ylabel('Kapur-Rokhlin error (bad quadrature)')
title(sprintf('$\\varepsilon = %1.1e$', epsi),'interpreter', 'latex')
axis tight

% end of file