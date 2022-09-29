%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% View the asymptotics of the vector potential integrand
% and test quadratures for the associated integral
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Set up shape and observation point
clear; clc; close all

version = 'normal'; % 'log' to see behavior near singularity

shape.kappa = 1;
shape.a = 1/3;
shape.R0 = 1;
shape.FB = 1;
shape.q0 = 1;
% Rmin = sqrt( shape.R0^2 - 2*shape.a*shape.R0 );
% Rmax = sqrt( shape.R0^2 + 2*shape.a*shape.R0 );

% t_zmax = roots([shape.a, shape.R0, shape.a]); % solve quadratic polynomial in cos(t_zmax)
% t_zmax = t_zmax(abs(t_zmax)<=1);
% t_zmax = acos(t_zmax);

% the observation point on the surface -- at t=t0
t0 = 1;
R = sqrt( shape.R0^2 + 2*shape.a*shape.R0*cos(t0) );
Z = shape.kappa*shape.a*shape.R0 ./ R .* sin(t0);

% KR convergence order parameter
order = 10;

% see how quadratures converge with number of quadrature nodes
nM = 40;
MM = round(exp(linspace(log(5), log(100), nM)))';

% table header
fprintf('Vector potential on surface (log-singular):\n')
fprintf('%6s | %20s | %20s\n', 'Nquad', 'Line (Alt. Trapz.)', 'Line (KR)')

% preallocate
[A_phi_altertrap, A_phi_KR] = deal(zeros(size(MM)));

%% Compare different quadratures for the singular integral

% main loop -- over number of quadrature points
for j = 1:nM
    
    M = MM(j);
    
    switch version
        
        case 'log' % for log-spaced parameter values near singularity
            M = max(MM);
            tt = logspace(-10, 0, 2*M)' + t0;
            [rr, zz, dpsi_dr, dpsi_dz, dr_dt, dz_dt, dpsi_dn_bracket, ...
                alpha, beta, k2] = get_functions_Jungpyo(shape, t0, tt);
            A_integrand = dpsi_dn_bracket .* 4 ./ sqrt(alpha + beta) ...
                .* ( 2./k2.*( ellipticK(k2) - ellipticE(k2) ) - ellipticK(k2) );
            break;

        otherwise % linear-spaced parameter values
            h = pi/M;
            tt = h.*(-M+1 : M)' + t0; %(-pi+h), ..., pi shifted by t0
    end
    
    % alternating trapezoidal
    if strcmpi(version,'normal')
        [rr, zz, dpsi_dr, dpsi_dz, dr_dt, dz_dt, dpsi_dn_bracket, ...
            alpha, beta, k2] = get_functions_Jungpyo(shape, t0, tt-h/2); % note the staggered grid "h/2"
        A_integrand = 1/(4*pi) .* dpsi_dn_bracket .* 4 ./ sqrt(alpha + beta) ...
            .* ( 2./k2.*( ellipticK(k2) - ellipticE(k2) ) - ellipticK(k2) );
        A_phi_altertrap(j) = h*sum(A_integrand);
    else
        A_phi_altertrap(j) = NaN;
    end

    % Kapur Rokhlin
    if M >= order+1 && strcmpi(version,'normal') % KR quadrature corrections are defined

        % generate integrand at equispaced points
        [rr, zz, dpsi_dr, dpsi_dz, dr_dt, dz_dt, dpsi_dn_bracket, ...
            alpha, beta, k2] = get_functions_Jungpyo(shape, t0, tt);
        A_integrand = 1/(4*pi) .* dpsi_dn_bracket .* 4 ./ sqrt(alpha + beta) ...
            .* ( 2./k2.*( ellipticK(k2) - ellipticE(k2) ) - ellipticK(k2) );

        % integrate by KR corrections
        A_phi_KR(j) = periodic_KR(h, A_integrand, order);

        % print to table
        fprintf( '%6d | %1.14e | %1.14e\n', ...
            2*M, A_phi_altertrap(j), A_phi_KR(j))
    else
        A_phi_KR(j) = NaN;
        fprintf( '%6d | %20s | %20s\n', ...
            2*M, 'N/A', 'N/A' )
    end
end

%% plot integrand to view asymptotics
figure(2); clf;
switch version
    case 'log'
        loglog(tt-t0, abs(A_integrand), '-')
        hold on;
        plot(tt-t0, abs(log(tt-t0)), 'k--')
        xlabel('$t - t_0$','interpreter','latex')
        legend('$\mathcal A$', 'log','interpreter','latex')
    otherwise
        plot(tt, A_integrand, '-')
        xlabel('$t$','interpreter','latex')
end
ylabel('$\mathcal A(t; \vec x)$','interpreter','latex')
title(['$A(\vec x) = \int_0^L \mathcal A(t; \vec x) \, dt$, ', sprintf('($t_0=%1.2f$)', t0) ], 'interpreter', 'latex')
axis tight
% saveas(gcf, sprintf('./img/vector_potential_integrand_t0_%d.png',t0))

%% self convergence test
if t0==1
    A_phi_star = -1.79169776399742e+00 / (4*pi);  % from KR with M=1e4, t0=1
elseif t0==0
    A_phi_star = -1.759217533099384e+00 / (4*pi); % ^ using t0=0
end
figure(1); clf;

relerr_altertrap = abs(A_phi_altertrap - A_phi_star)./abs(A_phi_star);
loglog(MM, relerr_altertrap, 'bx-')
hold on
relerr_KR = abs(A_phi_KR - A_phi_star)./abs(A_phi_star);
loglog(MM, relerr_KR, 'ro-')
xlabel('$2M$ (\# quadrature nodes)','interpreter','latex')
ylabel('Error (truth=KR w/ $M$=1e4)','interpreter','latex')

% asymptotics
loglog(MM, 1./MM, 'k--')
loglog(MM(MM>order), 2*max(relerr_KR).*(order./MM(MM>order)).^order, 'k:')
legend({'Alt. Trapz.', sprintf('%dth order KR',order), '1st order', sprintf('%dth order',order)})
title(sprintf('$t_0=%1.2f$', t0),'interpreter','latex')
ylim([min(relerr_KR), 1e1])
xlim([min(MM), max(MM)])
% saveas(gcf, sprintf('./img/vector_potential_error_convergence_t0_%d.png',t0))