%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Analyze the limiting behavior of a coefficient that appears when 
% reducing the harmonic function identity
% 1 + (1/2pi) * \int_\Gamma { dot(n(y), x-y) / || x-y ||^3 } dA(y) = 0
% (for x \in \Gamma) to a 1D integral
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; close all; clc

which_shape = 'Jungpyo';

%% set geometry

switch which_shape
    case 'Jungpyo'
        shape.kappa = 1;
        shape.a     = 1/3;
        shape.R0    = 1;
        shape.FB    = 1;
        shape.q0    = 1;
        figure; plot_surface_Jungpyo(shape, 100); hold on
        
        % computations if you want t = t_zmax, at the "top" of the generating curve
        t_zmax = roots([shape.a, shape.R0, shape.a]);   % solve quadratic polynomial in cos(t_zmax)
        t_zmax = t_zmax(abs(t_zmax)<=1);                % root which makes sense as a real cosine
        t_zmax = acos(t_zmax);

    case 'Andras'
        error('Andras geometry not yet supported.')
end

%% choose evaluation point (on surface)

% parameter of the evaluation point
t0 = 1;

switch which_shape
    case 'Jungpyo'
        [R, Z] = get_functions_Jungpyo(shape, t0, t0);
        hold on; plot(R,Z,'kx')
    case 'Andras'
        Rmax = 1+shape.epsi;
        R = Rmax;
        Z = 0;
end

%% set up discretized shape

M = 5e2; % sufficiently large number of parameter values
h = pi/M;
tt = h.*(-M+1 : M)' + t0;

switch which_shape
    case 'Jungpyo'
        [rr, zz, dpsi_dr, dpsi_dz, dr_dt, dz_dt, dpsi_dn_bracket, ...
            alpha, beta, k2, d2r_dt2, d2z_dt2] = get_functions_Jungpyo(shape, t0, tt);
%         d3r_dt3 = -shape.a*shape.R0 * ( -sin(tt) ./ rr - cos(tt) ./ rr.^2.*dr_dt ...
%             - cos(tt)./rr.^2.*dr_dt + 2*sin(tt)./rr.^3.*dr_dt.^2 - sin(tt)./rr.^2.*d2r_dt2 );
        
        % Fourier differentiation for numerical 3rd parameter deriv's
        d3r_dt3 = ifft(1i*[0:M-1, 0, -M+1:-1]'.*fft(d2r_dt2));
        d3z_dt3 = ifft(1i*[0:M-1, 0, -M+1:-1]'.*fft(d2z_dt2));

    case 'Andras'
        error('Andras geometry not yet supported.')
end

%% analysis -- limiting value of coefficient

coef = dz_dt .* (R - rr) - dr_dt .* (Z - zz);
coef = coef ./ ( alpha-beta );
coef_limit = 0.5 * (d2r_dt2(M) * dz_dt(M) - dr_dt(M) * d2z_dt2(M)) ./ (dr_dt(M)^2 + dz_dt(M)^2);

figure;
plot(tt, coef)
line([tt(1), tt(end)], [coef_limit, coef_limit])
line([t0, t0], [min(coef), max(coef)]) % crosshairs should intersect curve
% correct!

%% first derivative limit
coef1num = (((rr - R).*dz_dt - (zz - Z).*dr_dt).*(2*(rr - R).*dr_dt + 2*(zz - Z).*dz_dt)) - ((rr - R).*d2z_dt2 - (zz - Z).*d2r_dt2).*((rr - R).^2 + (zz - Z).^2);
coef1den = ((rr - R).^2 + (zz - Z).^2).^2;
coef1 = coef1num./coef1den;
coef1_limit = ( 2*(d2z_dt2(M)*dr_dt(M) - d2r_dt2(M)*dz_dt(M))*(6*dr_dt(M)*d2r_dt2(M) + 6*dz_dt(M)*d2z_dt2(M)) - 2*(2*dr_dt(M)^2 + 2*dz_dt(M)^2)*(2*dr_dt(M)*d3z_dt3(M) - 2*dz_dt(M)*d3r_dt3(M)) ) ...
    ./ ( 6*(2*dr_dt(M)^2 + 2*dz_dt(M)^2)^2 );

figure;
plot(tt, coef1)
line([tt(1), tt(end)], [coef1_limit, coef1_limit])
line([t0, t0], [min(coef1), max(coef1)]) % crosshairs should intersect curve
% correct!

%% 2nd derivative
coef2 = (2*((rr - R).*d2z_dt2 - (zz - Z).*d2r_dt2).*(2*(rr - R).*dr_dt + 2*(zz - Z).*dz_dt))./((rr - R).^2 + (zz - Z).^2).^2 ...
    - (d2z_dt2.*dr_dt - d2r_dt2.*dz_dt + (rr - R).*d3z_dt3 - (zz - Z).*d3r_dt3)./((rr - R).^2 + (zz - Z).^2) ...
    + (((rr - R).*dz_dt - (zz - Z).*dr_dt).*(2*(rr - R).*d2r_dt2 + 2*dr_dt.^2 + 2*dz_dt.^2 + 2*(zz - Z).*d2z_dt2))./((rr - R).^2 + (zz - Z).^2).^2 ...
    - (8*((rr - R).*dr_dt + (zz - Z).*dz_dt).^2.*((rr - R).*dz_dt - (zz - Z).*dr_dt))./((rr - R).^2 + (zz - Z).^2).^3;

figure;
plot(tt, coef2)