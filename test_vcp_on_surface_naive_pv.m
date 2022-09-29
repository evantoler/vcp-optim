%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Display asymptotics of a naive approach to the virtual casing principle
% on surface. Show that a principal value procedure is necessary.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; close all; clc

version = 'normal'; % 'log' to see behavior near boundary

% set shape parameters -- Jungpyo geometry
shape.kappa = 1;
shape.a = 1/3;
shape.R0 = 1;
shape.FB = 1;
shape.q0 = 1;
Rmin = sqrt( shape.R0^2 - 2*shape.a*shape.R0 );
Rmax = sqrt( shape.R0^2 + 2*shape.a*shape.R0 );

% calculations for picking a point at the "top" of the generating curve
t_crit = roots([shape.a, shape.R0, shape.a]); % solve quadratic polynomial in cos(t_zmax)
t_crit = t_crit(abs(t_crit)<=1);
t_crit = acos(t_crit);

% generate an observation point on surface at parameter t=t0
t0 = 1;
[R, Z, ~, ~, dr_dt0, dz_dt0, ~, ~, ~, ~] = get_functions_Jungpyo(shape, t0, t0);

% KR convergence order parameter
order = 10;

% see how quadratures converge with number of quadrature nodes
nM = 40;
MM = round(exp(linspace(log(5), log(500), nM)))';

% table header
fprintf('Vector potential on surface (log-singular):\n')
fprintf('%6s | %21s | %21s\n', 'Nquad', 'Line (Alt. Trapz.)', 'Line (KR)')

% preallocate
[B_normal_altertrap, B_normal_KR] = deal(zeros(size(MM)));

for j = 1:nM
    
    M = MM(j);
    
    %% Parameterize the surface
    
    switch version
        case 'log' % just for integrand -- ignore quadrature
            
            M = max(MM);
            tt = logspace(-10, 0, 2*M)' + t0;
            [rr, zz, dpsi_dr, dpsi_dz, dr_dt, dz_dt, dpsi_dn_bracket, ...
                alpha, beta, k2] = get_functions_Jungpyo(shape, t0, tt);
            
            % compute integrand
            B_n_integrand = -1/R * 2 / sqrt(dr_dt0^2 + dz_dt0^2) .* dpsi_dn_bracket .* 4 ./ sqrt(alpha + beta);
            term1 = -dr_dt0 * (-alpha./beta.*ellipticK(k2) + (alpha+beta)./beta.*ellipticE(k2));
            term2 = -R * (R*dr_dt0 + (Z-zz).*dz_dt0);
            term2 = term2 .* (-1./beta.*ellipticK(k2) + (1./(alpha-beta) + 1./beta).*ellipticE(k2));
            term3 = R*dr_dt0.*rr .* ( (2./beta - 1./(alpha-beta)).*ellipticE(k2) ...
                - 2./beta.*ellipticK(k2) + (alpha+beta)./(beta.^2).*(2-k2-2.*sqrt(1-k2)) );
            B_n_integrand = B_n_integrand .* (term1 + term2 + term3);
            
            break;
            
        otherwise
            h = pi/M;
            tt = h.*(-M+1 : M)' + t0; %(-pi+h), ..., pi shifted by t0
    end
    
    % alternating trapezoidal rule
    [rr, zz, dpsi_dr, dpsi_dz, dr_dt, dz_dt, dpsi_dn_bracket, ...
        alpha, beta, k2] = get_functions_Jungpyo(shape, t0, tt-h/2);
    
    B_n_integrand = -1/R * 2 / sqrt(dr_dt0^2 + dz_dt0^2) .* dpsi_dn_bracket .* 4 ./ sqrt(alpha + beta);
    term1 = -dr_dt0 * (-alpha./beta.*ellipticK(k2) + (alpha+beta)./beta.*ellipticE(k2));
    term2 = -R * (R*dr_dt0 + (Z-zz).*dz_dt0);
    term2 = term2 .* (-1./beta.*ellipticK(k2) + (1./(alpha-beta) + 1./beta).*ellipticE(k2));
    term3 = R*dr_dt0.*rr .* ( (2./beta - 1./(alpha-beta)).*ellipticE(k2) ...
        - 2./beta.*ellipticK(k2) + (alpha+beta)./(beta.^2).*(2-k2-2.*sqrt(1-k2)) );
    B_n_integrand = B_n_integrand .* (term1 + term2 + term3);
    
    B_normal_altertrap(j) = h*sum(B_n_integrand);    
    
    % Kapur Rokhlin rule
    if M >= order+1 % KR quadrature corrections are defined
        
        [rr, zz, dpsi_dr, dpsi_dz, dr_dt, dz_dt, dpsi_dn_bracket, ...
            alpha, beta, k2] = get_functions_Jungpyo(shape, t0, tt);
        
        B_n_integrand = -1/R * 2 / sqrt(dr_dt0^2 + dz_dt0^2) .* dpsi_dn_bracket .* 4 ./ sqrt(alpha + beta);
        term1 = -dr_dt0 * (-alpha./beta.*ellipticK(k2) + (alpha+beta)./beta.*ellipticE(k2));
        term2 = -R * (R*dr_dt0 + (Z-zz).*dz_dt0);
        term2 = term2 .* (-1./beta.*ellipticK(k2) + (1./(alpha-beta) + 1./beta).*ellipticE(k2));
        term3 = R*dr_dt0.*rr .* ( (2./beta - 1./(alpha-beta)).*ellipticE(k2) ...
            - 2./beta.*ellipticK(k2) + (alpha+beta)./(beta.^2).*(2-k2-2.*sqrt(1-k2)) );
        B_n_integrand = B_n_integrand .* (term1 + term2 + term3);
        
        B_normal_KR(j) = periodic_KR(h, B_n_integrand, order);

        % print table line
        fprintf( '%6d | %+1.14e | %1.14e\n', ...
            2*M, B_normal_altertrap(j), B_normal_KR(j))
        
    else
        B_normal_KR(j) = NaN;
        fprintf( '%6d | %+1.14e | %20s\n', ...
            2*M, B_normal_altertrap(j), 'N/A' )
    end
end

%% plot integrand
figure;
switch version
    case 'log'
        loglog(tt-t0, abs(B_n_integrand), '-')
        hold on;
        plot(tt-t0, abs(1./(tt-t0).^2), 'k--')
        xlabel('$t - t_0$','interpreter','latex')
        legend('Integrand', '$1/|t-t_0|^2$','interpreter','latex')
    otherwise
        plot(tt, B_n_integrand, '-')
        xlabel('$t$','interpreter','latex')
end
ylabel 'Integrand'
title('Integrand of $\vec n \cdot \vec B$','interpreter','latex')
xlabel('$t$','interpreter','latex')
axis tight
% saveas(gcf, sprintf('./img/tangential_deriv_integrand_t0_%d.png',t0))