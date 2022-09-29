%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% evaluate magnetic field on-surface with Zakharov line integral
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all; clear; clc

shape.kappa = 1;
shape.a = 1/3;
shape.R0 = 1;
shape.FB = 1;
shape.q0 = 1;
Rmin = sqrt( shape.R0^2 - 2*shape.a*shape.R0 );
Rmax = sqrt( shape.R0^2 + 2*shape.a*shape.R0 );

t0 = 0;
R = sqrt( shape.R0^2 + 2*shape.a*shape.R0*cos(t0) );
Z = shape.kappa*shape.a*shape.R0 ./ R .* sin(t0);

%% visualize surface

% N = 1e3;
% figure(1); clf;
% plot_surface_Jungpyo(shape, N)
% hold on
% plot(R,Z,'kx')

%% compute "true" field

[B_rSTAR, B_zSTAR] = volume_integral_Jungpyo(R, Z, shape, 'builtin');

% alternative options: high-res Kapur-Rokhlin or line integral
% [B_rSTAR, B_zSTAR] = volume_integral_Jungpyo(R, Z, shape, 'KR', 1e3);
% [B_rSTAR, B_zSTAR] = Zakharov_line_integral_Jungpyo(r, z, shape, 1e3);

%% field as the target point approaches the boundary
% Uncomment portions to create a video

Ndr = 20; % number of positive "dr" perturbations to test

% v = VideoWriter('Line_Integrand_approach_surface.avi');
% v.FrameRate=Ndr;
% open(v);
% figure(1); clf;

% preallocate
[B_rV,B_zV,B_rLoff,B_zLoff,B_rLon,B_zLon,B_rS,B_zS] = deal( zeros(Ndr+1,1) );

M = 100; % fixed number of quadrature points
dr = [0, exp(linspace(log(1e-10), log(1e0), Ndr))]; % perturb the point off the boundary

% main loop: get closer to the surface
for j = length(dr):-1:1
    
    fprintf('Iter %d\n', j)
    
    [B_rV(j),   B_zV(j)]        = volume_integral_Jungpyo(R+dr(j), Z, shape, 'builtin');
    [B_rLoff(j),B_zLoff(j)]     = Zakharov_line_integral_Jungpyo(R+dr(j), Z, shape, M, false); % treat as off surface
    [B_rLon(j), B_zLon(j)]      = Zakharov_line_integral_Jungpyo(R+dr(j), Z, shape, M, true); % treat as on surface
    [B_rS(j),   B_zS(j)]        = surface_integral_Jungpyo(R+dr(j), Z, shape, M/2); % surface integral
    
%     title(sprintf('$B_z = 2 \\cdot \\frac 1{4\\pi} \\int_{-\\pi}^\\pi f(t) \\, dt$;   $(R, Z)=(R_{max} + %1.0e, 0)$', dr(j)))
%     writeVideo(v, getframe(gcf)); % using global v = VideoWriter
end

% close(v)

%% plotting

figure; clf;
subplot(1,2,1)
loglog(dr, abs(B_rV), 'bx:', dr, abs(B_rLoff), 'ro:', dr, abs(B_rLon), 'gs:')
hold on
loglog(dr, abs(B_rS), 'm*')
xlabel('$\Delta R$', 'interpreter','latex')
ylabel 'Field value'
title('$B_r(R+\Delta R, Z=0)$','interpreter','latex')
legend('Volume Int','Line Int (off $\Gamma$; no jump)','Line Int (on $\Gamma$; w/ jump)','Surface Int (on $\Gamma$; w/jump)','interpreter','latex')
subplot(1,2,2)
loglog(dr, abs(B_zV), 'bx-', dr, abs(B_zLoff), 'ro-', dr, abs(B_zLon), 'gs-')
hold on
loglog(dr, abs(B_zS), 'm*')
xlabel('$\Delta R$', 'interpreter','latex')
ylabel 'Field value'
title('$B_z(R+\Delta R, Z=0)$','interpreter','latex')
legend('Volume Int','Line Int (off $\Gamma$; no jump)','Line Int (on $\Gamma$; w/ jump)','Surface Int (on $\Gamma$; w/jump)','interpreter','latex')
% saveas(gcf, 'approach_Alt_on.fig');
% saveas(gcf, 'approach_Alt_on.png');

figure;
subplot(1,2,1)
loglog(dr, abs(B_rV - B_rSTAR), 'bx:', dr, abs(B_rLoff - B_rSTAR), 'ro:', dr, abs(B_rLon - B_rSTAR), 'gs:')
hold on
loglog(dr, abs(B_rS - B_rSTAR), 'm*')
xlabel('$\Delta R$', 'interpreter','latex')
ylabel 'Absolute Error'
title('Error in $B_r(R+\Delta R, Z=0)$','interpreter','latex')
legend('Volume Int','Line Int (off $\Gamma$; no jump)','Line Int (on $\Gamma$; w/ jump)','Surface Int (on $\Gamma$; w/jump)','interpreter','latex')
subplot(1,2,2)
loglog(dr, abs(B_zV - B_zSTAR)./abs(B_zSTAR), 'bx-', dr, abs(B_zLoff - B_zSTAR)./abs(B_zSTAR), 'ro-', dr, abs(B_zLon - B_zSTAR)./abs(B_zSTAR), 'gs-')
hold on
loglog(dr, abs(B_zS - B_zSTAR), 'm*')
xlabel('$\Delta R$', 'interpreter','latex')
ylabel 'Relative Error'
title('Error in $B_z(R+\Delta R, Z=0)$','interpreter','latex')
legend('Volume Int','Line Int (off $\Gamma$; no jump)','Line Int (on $\Gamma$; w/ jump)','Surface Int (on $\Gamma$; w/jump)','interpreter','latex')

%%

% precomputed figures
% open('approach_KR_on.fig');
% open('approach_Alt_on.fig');

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% convergence and time scaling with N
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% normalization constant for relative error
norm_Bstar = sqrt(B_rSTAR^2 + B_zSTAR^2);

order = 10; % Kapur-Rokhlin convergence order parameter
NN = round(linspace(order+1, 100/2, Ndr)); % order+1 necessary for KR
NN = 2*NN; % use even N for convenience
[errV, timeV, errZ, timeZ] ...
    = deal( zeros(size(NN)) ); % V=volume, Z=Zakharov

for j = length(NN):-1:1
    
    fprintf('Iter%d\n',j)
    N = NN(j);
    
    % perform volume integration
    tic;
    [B_rV, B_zV] = volume_integral_Jungpyo(R, Z, shape, 'gauss', N);
    timeV(j) = toc;
    errV(j) = norm([B_rV; B_zV] - [B_rSTAR; B_zSTAR])/norm_Bstar;
    
    % perform line integration
    tic;
    [B_rZ, B_zZ] = Zakharov_line_integral_Jungpyo(R, Z, shape, N/2, true, t0, false, order);
    timeZ(j) = toc;
    errZ(j) = norm([B_rZ; B_zZ] - [B_rSTAR; B_zSTAR])/norm_Bstar;
    
end

figure;
loglog(NN, errV, 'bx-', NN, errZ, 'ro-'); hold on
loglog(NN, 10*NN.^(-3), 'k--')
loglog(NN(NN<50), errZ(1)*(NN(1)./NN(NN<50)).^order, 'k:')
xlabel('$N$ (\# quadrature points)','fontsize',20,'interpreter','latex')
ylabel('$\|\vec B_{approx} - \vec B_*\|_2/\|\vec B_*\|_2$','fontsize',20,'interpreter','latex')
title('Jungpyo''s geometry','fontsize',20)
axis tight
legend('volume', 'Zakharov','$O(N^{-3})$', '$O(e^{-0.7N})$', ...
    'interpreter','latex', 'location', 'east')
% saveas(gcf, 'error_convergence_Zakharov_as_trapezoid_Jungpyo.jpg')

figure;
semilogy(NN, timeV, 'bx-', NN, timeZ, 'ro-')
xlabel('$N$ (\# quadrature points)','fontsize',20,'interpreter','latex')
ylabel('time (sec)','fontsize',20)
title('Jungpyo''s geometry','fontsize',20)
axis tight
legend('volume', 'Zakharov', 'location', 'northwest')

