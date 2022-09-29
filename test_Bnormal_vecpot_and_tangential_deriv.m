%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% examine the integrands involved with a weakly singular integral for the
% normal component of the magnetic field. We use the intermediary 
% vector potential and use the Zakharov line integral method.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; clc; close all

%% Set up shape and observation points

C = 1/(4*pi);

% set shape parameters
shape.kappa = 1;
shape.a = 1/3;
shape.R0 = 1;
shape.FB = 1;
shape.q0 = 1;
% Rmin = sqrt( shape.R0^2 - 2*shape.a*shape.R0 );
% Rmax = sqrt( shape.R0^2 + 2*shape.a*shape.R0 );

% KR convergence order parameter
order = 10;

% outer loop: over evaluation points x on the surface
% loop over t0 in [0, 2*pi]
Mx = 50;
hx = pi/Mx;
t0t0 = hx.*(-Mx+1 : Mx)';
[RR, ZZ, ~, ~, dR_dt0, dZ_dt0, ~, ~, ~, ~] = get_functions_Jungpyo(shape, 0, t0t0);

% inner discretization for generating curve
% fixed for each evaluation point
My = 100;
hy = pi/My;

% preallocate
[A_phi_altertrap, A_phi_KR] = deal(zeros(size(t0t0)));

%% Compute the vector potential 

for j = 1:2*Mx
    
    fprintf('Iter %d/%d\n', j, 2*Mx)
    t0 = t0t0(j);
    
    % Parameterize the surface
    
    tt = hy.*(-My+1 : My)' + t0; %(-pi+h), ..., pi shifted by t0
    
    % alternating trapezoidal rule
    [rr, zz, dpsi_dr, dpsi_dz, dr_dt, dz_dt, dpsi_dn_bracket, ...
        alpha, beta, k2] = get_functions_Jungpyo(shape, t0, tt-hy/2);
    A_integrand = C * dpsi_dn_bracket .* 4 ./ sqrt(alpha + beta) ...
        .* ( 2./k2.*( ellipticK(k2) - ellipticE(k2) ) - ellipticK(k2) );
    A_phi_altertrap(j) = hy*sum(A_integrand);
    
    % Kapur Rokhlin rule
    if My >= order+1 % KR quadrature corrections are defined
        [rr, zz, dpsi_dr, dpsi_dz, dr_dt, dz_dt, dpsi_dn_bracket, ...
            alpha, beta, k2] = get_functions_Jungpyo(shape, t0, tt);
        A_integrand = C * dpsi_dn_bracket .* 4 ./ sqrt(alpha + beta) ...
            .* ( 2./k2.*( ellipticK(k2) - ellipticE(k2) ) - ellipticK(k2) );
        A_phi_KR(j) = periodic_KR(hy, A_integrand, order);
    else
        A_phi_KR(j) = NaN;
    end
end

%% Plot integrand

figure;
plot(tt, A_integrand, '-')
xlabel('$t$','interpreter','latex')
ylabel('$\mathcal A(t; \vec x)$','interpreter','latex')
title(['$A(\vec x) = \int_0^L \mathcal A(t; \vec x) \, dt$, ', sprintf('($t_0=%1.2f$)', t0) ], 'interpreter', 'latex')
axis tight
% saveas(gcf, sprintf('./img/vector_potential_integrand_t0_%d.png',t0))

%% Plot vector potential on surface

figure;
plot(t0t0, A_phi_altertrap, 'b-')
hold on
plot(t0t0, A_phi_KR, 'r-')
xlabel('$t_0$','interpreter','latex')
ylabel('$\vec A(\vec y(t=t_0, \theta=0))$','interpreter','latex')
title(sprintf('$M_x = %d, M_y = %d$', Mx, My),'interpreter','latex')
legend('Alter. Trapz.','10th order KR')
axis tight

%% Take tangential derivative to get dot(n, B)

% "virtual casing" flux function for magnetic field just due to plasma
psi_VC = RR .* A_phi_KR;

% Fourier differentiation in t0
B_normal = ifft(1i*[0:Mx-1, 0, -Mx+1:-1]'.*fft(psi_VC));
B_normal = B_normal .* 2 ./ sqrt(dR_dt0.^2 + dZ_dt0.^2);    % chain rule factor
B_normal = -B_normal ./ RR;                                 % factor -1/R from triple product sign change and grad(e_phi)=(e_phi)/R

%% plot dot(n, B)

figure;
plot(t0t0, B_normal, 'x-')
xlabel('$t_0$','interpreter','latex')
ylabel('$\vec n \cdot \vec B$','interpreter','latex')
axis tight
