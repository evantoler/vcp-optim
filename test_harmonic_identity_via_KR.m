%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Test the following identity with Kapur-Rokhlin quadrature:
% 1 + (1/2pi) * \int_\Gamma { dot(n(y), x-y) / || x-y ||^3 } dA(y) = 0
% for x \in \Gamma
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; close all; clc

which_shape = 'Jungpyo'; % alternatively, 'Andras'

%% set geometry

switch which_shape
    case 'Jungpyo'
        shape.kappa = 1;
        shape.a     = 1/3;
        shape.R0    = 1;
        shape.FB    = 1;
        shape.q0    = 1;
        figure; plot_surface_Jungpyo(shape, 100); hold on
        
    case 'Andras'
        shape.kappa = 2;
        shape.delta = 0.7;
        shape.epsi = 0.2;
        figure; plot_surface_Andras(shape, 100); hold on
end

%% choose evaluation point (on surface)

% computations if you want t = t_zmax, at the "top" of the generating curve
t_zmax = roots([shape.a, shape.R0, shape.a]);   % solve quadratic polynomial in cos(t_zmax)
t_zmax = t_zmax(abs(t_zmax)<=1);                % root which makes sense as a real cosine
t_zmax = acos(t_zmax);

t0 = 1;

switch which_shape
    case 'Jungpyo'
        [R, Z] = get_functions_Jungpyo(shape, t0, t0);
    case 'Andras'
        Rmax = 1+shape.epsi;
        R = Rmax;
        Z = 0;
end

hold on; plot(R,Z,'kx')

%% set up convergence test

num_M = 30;
order = 10; % KR convergence order parameter

% using 2*M quadrature nodes
MM = round( exp(linspace(log(order+1), log(2e2), num_M)) )';

% preallocate
% {punctured, KR, alternating trapezoidal, surface} quadrature rules
[err, errKR, errAlt, err2D] = deal( zeros(size(MM)) );

%% main convergence loop

wb = waitbar(0, 'Computing double layer test...');
for j = 1:num_M
    
    %% discretize points on boundary
    
    % 2*M equispaced parameters
    M = MM(j);
    h = pi/M;
    tt = h.*(-M+1 : M)' + t0;
    
    switch which_shape
        case 'Jungpyo'
            [rr, zz, ~, ~, dr_dt, dz_dt, ~, alpha, beta, k2, d2r_dt2, d2z_dt2] ...
                = get_functions_Jungpyo(shape, t0, tt);
        case 'Andras'
            error('Andras geometry not yet supported.')
    end
    
    %% compute 1D integrand (after analytically integrating out angle)
    
    coef = ( dz_dt.*(R-rr) - dr_dt.*(Z-zz) ) ./ (1-k2);
    ff = -2*R .* dz_dt ./ k2 .* ellipticK(k2) ...
        + ( 2*R .* dz_dt ./ k2 + coef ) .* ellipticE(k2);
    ff = ff .* 4 .* rr ./ (alpha + beta).^(3/2);
    
    %% compute 2D surface integrand
    [T, THETA] = meshgrid(tt, h.*(-M+1 : M)' );
    [RR, ZZ, ~, ~, DR_DT, DZ_DT, ~, ALPHA, BETA] = get_functions_Jungpyo(shape, t0, T);
    f2D = DZ_DT.* R .* cos(THETA) - ( DZ_DT.*RR + DR_DT.*(Z - ZZ) );
    f2D = f2D .* RR ./ (ALPHA - BETA .* cos(THETA)).^(3/2);
    f2D(M,M) = 0; % not defined at x -- omit
    
    %% integrate; apply KR correction
    
    % punctured trapezoidal rule
    I = h*sum(ff([1:M-1, M+1:end]));

    % KR correction
    IKR = periodic_KR(h,ff,order);
    
    % alternating trapezoidal
    IAltertrap = (2*h) * sum( ff([M-1:-2:1, M+1:2:end]) );
    
    % 2D
    I2D = h^2 * sum(f2D, 'all');
    
    %% track error
    
    I = 1 + I/(2*pi); IKR = 1 + IKR/(2*pi); IAltertrap = 1 + IAltertrap/(2*pi); I2D = 1 + I2D/(2*pi);
    err(j) = abs(I - 0); errKR(j) = abs(IKR - 0); errAlt(j) = abs(IAltertrap - 0); err2D(j) = abs(I2D - 0);
    
    waitbar(j/num_M, wb)
end

close(wb)

%% debugging -- check limiting value of coef
% coef_limit = 2*R^2 * (d2r_dt2(M) * dz_dt(M) - dr_dt(M) * d2z_dt2(M)) ./ (dr_dt(M)^2 + dz_dt(M)^2);
% figure;
% plot(tt, coef)
% line([tt(1), tt(end)], [coef_limit, coef_limit])
% line([t0, t0], [min(coef), max(coef)]) % crosshairs should intersect curve
% % correct!

%% plot integrand
% 
% figure;
% plot(tt, ff)
% xlabel('$t$','interpreter','latex')
% ylabel('1D Integrand')
% axis tight

%% plot convergence

figure;
loglog(2*MM-1, errKR, 'bx-', 2*MM-1, errAlt, 'ro--')
% hold on; loglog((2*MM).^2-1, err2D, 'g*') % uncomment to include surface integral
xlabel 'Number of quadrature nodes'
ylabel 'Quadrature error'
title(sprintf('Harmonic Function Test at $t_0 = %1.1f$', t0),'interpreter','latex')
hold on
loglog(2*MM(MM<100)-1, errKR(1)*(MM(1)./MM(MM<100)).^(order), 'k:')
legend('Kapur-Rokhlin', 'Alternating Trapezoidal', sprintf('$O(h^{%d})$',order),'location','southwest', 'interpreter', 'latex')
axis tight
ax = gca;
ax.XAxis.TickValues = sort(unique( [2*MM(1)-1, 2*MM(end)-1, logspace(0, 6, 7)] ));
% saveas(gcf, sprintf('./img/harmonic_test_t0_%d.eps', t0), 'epsc')