function plot_surface_Jungpyo(shape, N)
% shape   (struct) geometry parameters with fields as outlined 
%             in Lee & Cerfon:
%             .kappa
%             .a
%             .R0
%             .FB
%             .q0
% N       (int >0) number of points to use in discretizing hte surface

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot the surface
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% parameterize the generating curve
h = 2*pi/N;
t = h.*(1:N)';
r = sqrt( shape.R0^2 + 2*shape.a*shape.R0*cos(t) );
z = shape.kappa*shape.a*shape.R0 ./ r .* sin(t);

% visualize surface
plot(r, z, 'b-')
axis equal
xlabel $r$
ylabel $z$

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% visually check that bounds are correct
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% hold on
% Rmin = sqrt( shape.R0^2 - 2*shape.a*shape.R0 );
% Rmax = sqrt( shape.R0^2 + 2*shape.a*shape.R0 );
% rr   = rprime(1:ceil(N/2)); % get rprime values from the upper half
% Zmax = shape.kappa*shape.a*shape.R0./rr .* sqrt(1 - ((rr.^2-shape.R0^2)./(2*shape.a*shape.R0)).^2);
% % Zmax = real(Zmax);
% Zmin = -Zmax; % Z bounds
% plot(rr, Zmin, 'r:', rr, Zmax, 'g:')
% line([Rmin Rmin], [min(Zmin) max(Zmax)])
% line([Rmax Rmax], [min(Zmin) max(Zmax)])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% check the normal vector direction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% hold on
% drprime_dt = -shape.a*shape.R0 ./ rprime .* sin(t);
% dzprime_dt = shape.kappa*shape.a*shape.R0 ...
%     .* ( cos(t)./rprime + shape.a*shape.R0.*(sin(t)).^2./(rprime.^3) );
% 
% Jac =  rprime .* sqrt(drprime_dt.^2 + dzprime_dt.^2);
% n_r =  rprime ./ Jac .* dzprime_dt;
% n_z = -rprime ./ Jac .* drprime_dt;
% quiver(rprime, zprime, n_r, n_z)

return