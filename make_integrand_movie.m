%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% examine the integrands involved with the virtual casing principle
% applied on surface. We use the Zakharov line integral method.
% This creates an integral requiring principal value analysis.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Uncomment sections to export video file

clear; close all; clc

%% Set up shape and observation points

M = 100;            % number of discretization points = 2*M
epsi = 7.5e-2;      % p.v. parameter
version = 'normal'; % 'log' to see behavior near boundary
nt0 = 200;

% set shape parameters -- Jungpyo geometry
shape.kappa = 1;
shape.a = 1/3;
shape.R0 = 1;
shape.FB = 1;
shape.q0 = 1;
Rmin = sqrt( shape.R0^2 - 2*shape.a*shape.R0 );
Rmax = sqrt( shape.R0^2 + 2*shape.a*shape.R0 );

omega = (pi-epsi)/2;
t0t0 = 2*pi * (0:nt0-1)'./nt0;

% make sure to include points with maximum/minimum Z coordinate
t_zmax = roots([shape.a, shape.R0, shape.a]); % solve quadratic polynomial in cos(t_zmax)
t_zmax = t_zmax(abs(t_zmax)<=1);
t_zmax = acos(t_zmax);
t0t0 = [t0t0; t_zmax; 2*pi-t_zmax]; % one critical t value maxes Z; the other, min
t0t0 = sort(unique(t0t0));

% determine Zmax and Zmin
Rcrit = sqrt( shape.R0^2 + 2*shape.a*shape.R0*cos(t_zmax) );
Zmax = shape.kappa*shape.a*shape.R0 ./ Rcrit .* sin(t_zmax);
Rcrit = sqrt( shape.R0^2 + 2*shape.a*shape.R0*cos(2*pi-t_zmax) );
Zmin = shape.kappa*shape.a*shape.R0 ./ Rcrit .* sin(2*pi-t_zmax);
clear Rcrit

fig = figure(1);
fig.Position(3:4) = [800   600]; % make 600px wide & 800 px tall

% v = VideoWriter('./img/integrand_movie.avi');
% v.FrameRate=5;
% open(v);

C = 1/(4*pi); % common constant factor

for j = 1 : length(t0t0)
    
    %% Parameterize surface
    
    % compute observation point
    t0 = t0t0(j);
    R = sqrt( shape.R0^2 + 2*shape.a*shape.R0*cos(t0) );
    Z = shape.kappa*shape.a*shape.R0 ./ R .* sin(t0);
    
    % make discretized parameters
    switch version
        case 'log'
            tt = logspace(-10, 0, 2*M)' + t0;
        otherwise
            h = pi/M;
            tt = h.*(-M+1 : M)' + t0; %(-pi+h), ..., pi shifted by t0
    end
    
    % compute important quantities
    [rr, zz, dpsi_dr, dpsi_dz, dr_dt, dz_dt, dpsi_dn_bracket, ...
        alpha, beta, k2, d2r_dt2, d2z_dt2] = get_functions_Jungpyo(shape, t0, tt);
    
    %% define integrands in respective regions
    
    % determine which points are near the singularity
    idx_small = ( abs(tt - t0) <= epsi );
    
    % unadjusted
    % Iz = C * dpsi_dn_bracket * 4./(alpha+beta).^(3/2);
    
    % for t away from singularity at t0
    fz_far = C * dpsi_dn_bracket * 4./(alpha+beta).^(3/2) ...
        .* ( ...
        ( (rr - R)./(1-k2) - 2*R./k2 ) .* ellipticE(k2) ...
        + 2*R./k2 .* ellipticK(k2) ...
        );
    fz_far(idx_small) = 0;
    
    fr_far = C * dpsi_dn_bracket * 4./(alpha+beta).^(3/2) ...
        .* (Z - zz) ...
        .* ( ...
        ( 2./k2 + 1./(1-k2) ) .* ellipticE(k2) ...
        - 2./k2 .* ellipticK(k2) ...
        );
    fr_far(idx_small) = 0;
    
    % debugging -- check limiting values
    % for t near the singularity at t0
%     fz_near = C * dpsi_dn_bracket * 4./(alpha+beta).^(3/2) ...
%         .* ( ...
%         ( (rr - R)./(1-k2) - 2*R./k2 ) .* ellipticE(omega, k2) ...
%         + 2*R./k2 .* ellipticF(omega, k2) ...
%         + ( (R-rr).*k2./(1-k2) + 2*R ) ...
%         .* sin(omega)*cos(omega)./sqrt(1-k2.*(sin(omega)^2)) ...
%         );
%     fz_near(~idx_small) = 0;
%     fz_near(M) = C * dpsi_dn_bracket(M) * 4./(alpha(M)+beta(M)).^(3/2) ...
%         .* ( ...
%         ( 2*R^2*d2r_dt2(M)/dz_dt(M)^2 - 2*R./k2(M) ) .* ellipticE(omega, k2(M)) ...
%         + 2*R./k2(M) .* ellipticF(omega, k2(M)) ...
%         + ( -2*R^2*d2r_dt2(M)/dz_dt(M)^2 .* k2(M) + 2*R ) ...
%         .* sin(omega)*cos(omega)./sqrt(1-k2(M).*(sin(omega)^2)) ...
%         ); % limiting value as t -> t0; removable singularity
    
    %% Plot integrands into movie
    
    figure(1); clf
    
    subplot(1,2,1)
    
    plot_surface_Jungpyo(shape, 1000);
    hold on;
    plot(R, Z, 'kx')
    xlim([Rmin-0.01, Rmax + 0.01])
    ylim(1.01*[Zmin, Zmax])
    title(sprintf('$t_0 = %1.2f$', t0),'interpreter', 'latex')
    
    subplot(1,2,2)
    
    switch version
        case 'log'
            legendstr = {'$|t-t_0|>\varepsilon$', '$|t-t_0|\le \varepsilon$', '$1/|t-t_0|$'};
            loglog(tt-t0, abs(fz_far), tt-t0, abs(fz_near), tt-t0, abs(1./(tt-t0)),'k--')
            xlabel('$t - t_0$','interpreter','latex')
        otherwise
            legendstr = {'$f_r$', '$f_z$'};
            plot(tt, fr_far)
            hold on
            plot(tt, fz_far)
            % debugging for limiting part
%             hold on;
%             plot(t, fz_near)
%             legendstr = {'$|t-t_0|>\varepsilon$'}; % , '$|t-t_0|\le \varepsilon$'};
            xlabel('$t$','interpreter','latex')
            xlim([-pi, pi] + t0)
            ylim([-0.8, 0.8]) % hard coded to keep consistent across frames
    end
    
    ylabel('$f(t; \epsilon)$ components','interpreter', 'latex')
    title(sprintf('$\\varepsilon = %1.1e$', epsi),'interpreter', 'latex')
    legend(legendstr, 'interpreter', 'latex', 'location', 'northeast')
    
    % video
%     writeVideo(v, getframe(gcf));
    drawnow
    
end

% close(v)