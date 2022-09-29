%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Compare different quadrature rules for computing the magnetic field
%  on a plasma flux surface, caused by circulating plasma
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; clc

%% Set up shape and observation point

% Jungpyo geometry parameters
shape.kappa = 1;
shape.a = 1/3;
shape.R0 = 1;
shape.FB = 1;
shape.q0 = 1;
Rmin = sqrt( shape.R0^2 - 2*shape.a*shape.R0 );
Rmax = sqrt( shape.R0^2 + 2*shape.a*shape.R0 );

% generate observation point on surface
t0 = 1;
fprintf('t0=%1.1f\n', t0)
R = sqrt( shape.R0^2 + 2*shape.a*shape.R0*cos(t0) );
Z = shape.kappa*shape.a*shape.R0 ./ R .* sin(t0);

% Kapur-Rokhlin parameters
order = 10;

% print table of convergence --- header line
fstr = '%1.15f';
fprintf('B_z on surface:\n')
fprintf('%6s | %20s | %20s | %20s | %20s\n', 'Nquad', 'Volume', 'Surface', 'Line (Alt. Trapz.)', 'Line (KR)')

%% Compute field for different numbers of quadrature points

nM = 40;
MM = round(exp(linspace(log(5), log(500), nM)))';

[B_rV,B_zV,B_rLAlt,B_zLAlt,B_rLKR,B_zLKR,B_rS,B_zS] = deal( zeros(size(MM)) );

for j = 1:nM
    
    M = MM(j);
    
    % volume integral using Gauss quadrature
    [B_rV(j), B_zV(j)] ...
        = volume_integral_Jungpyo(R, Z, shape, 'gauss', 2*M);

    % line integral on surface at t0 using alternating trapezoidal rule
    [B_rLAlt(j), B_zLAlt(j)]    ...
        = Zakharov_line_integral_Jungpyo(R, Z, shape, M, true, t0, true);
    
    % surface integral using 2D alternating trapezoidal rule
    [B_rS(j), B_zS(j)] ...
        = surface_integral_Jungpyo(R, Z, shape, M, true, t0, false);
    
    % Kapur-Rokhlin rule on the 1D reduction of the virtual casing integral
    if M >= order+1 % KR quadrature corrections are defined
        [B_rLKR(j), B_zLKR(j)]  ...
            = Zakharov_line_integral_Jungpyo(R, Z, shape, M, true, t0, false, order);

        % print table entry
        fprintf( '%6d | %1.14e | %1.14e | %1.14e | %1.14e\n', ...
            2*M, B_zV(j), B_zS(j), B_zLAlt(j), B_zLKR(j) )
    else
        fprintf( '%6d | %1.14e | %1.14e | %1.14e | %20s\n', ...
            2*M, B_zV(j), B_zS(j), B_zLAlt(j), 'N/A' )
    end
    
end
