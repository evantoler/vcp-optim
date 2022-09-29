%% 0. Set up shape and observation point
clear; clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% examine the integrands involved with a principal value integral
% using the Zakharov line integral method
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

M = 1000;               % number of discretization points = 2*M
% epsi = 1e-2;          % p.v. parameter
version = 'normal';     % 'log' to see behavior near boundary

% set shape parameters
shape.kappa = 1;
shape.a = 1/3;
shape.R0 = 1;
shape.FB = 1;
shape.q0 = 1;
Rmin = sqrt( shape.R0^2 - 2*shape.a*shape.R0 );
Rmax = sqrt( shape.R0^2 + 2*shape.a*shape.R0 );

% find parameter value corresponding to maximum value of z
t_crit = roots([shape.a, shape.R0, shape.a]);   % solve quadratic polynomial in cos(t_zmax)
t_crit = t_crit(abs(t_crit)<=1);
t_crit = acos(t_crit);

% set observation point
t0 = t_crit;
R = sqrt( shape.R0^2 + 2*shape.a*shape.R0*cos(t0) );
Z = shape.kappa*shape.a*shape.R0 ./ R .* sin(t0);

% omega = (pi-epsi)/2;

%% 1. Surface; parametrize
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% generate surface and important quantities
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch version
    case 'log'
        t = logspace(-10, 0, 2*M)' + t0;
    otherwise
        h = pi/M;
        t = h.*(-M+1 : M)' + t0; %(-pi+h), ..., pi shifted by t0
end

rr = sqrt( shape.R0^2 + 2*shape.a*shape.R0*cos(t) );
zz = shape.kappa*shape.a*shape.R0 ./ rr .* sin(t);

% dpsi_dr = shape.kappa*shape.FB/(2*shape.R0^3*shape.q0) .* rr ...
%     .* ( rr.^2-shape.R0^2 + 2*zz.^2./shape.kappa^2 );
% dpsi_dz = shape.a*shape.FB/(shape.R0^2*shape.q0) .* rr .* sin(t);
% 
% dr_dt = -shape.a*shape.R0*sin(t) ./ rr;
% dz_dt = shape.kappa*shape.a*shape.R0 * (cos(t)./rr - sin(t)./rr.^2.*dr_dt);
% 
% d2r_dt2 = -shape.a*shape.R0* (cos(t) ./ rr - sin(t)./rr.^2.*dr_dt);
% 
% dpsi_dn_bracket = dpsi_dz.*dr_dt - dpsi_dr.*dz_dt;

% see: ellipticE is uniformly continuous near t0
alpha = R^2 + rr.^2 + (Z-zz).^2;
beta = 2*R.*rr;
k2 = 2*beta ./ ( alpha + beta );
figure; plot(t, ellipticE(k2), 'b-')

