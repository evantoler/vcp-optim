%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Verify that the vector potential approach to computing dot(n,B) on the
% surface agrees with results from directly using the Biot-Savart Law.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; clc; close all

%% Parameterize surface

% shape parameters
shape.kappa = 1;
shape.a = 1/3;
shape.R0 = 1;
shape.FB = 1;
shape.q0 = 1;

% discretization numbers
Mx = 20;

% observation point discretization
hx = pi/Mx;
t0t0 = hx.*(-Mx : Mx-1)' + pi; % [-pi, pi] + pi
[RR, ZZ, ~, ~, dR_dt0, dZ_dt0, ~, ~, ~, ~] = get_functions_Jungpyo(shape, 0, t0t0); % dummy value of t0=0

% KR parameter
order = 10;

%% Method 1: Use vector potential

My = 50;
B_normal_vecpot = B_normal_Jungpyo(shape, t0t0, My, order);
fprintf('Vecpot done.\n')

%% Methods 2,3,4,5: "normal dotted with full vector result" method

% preallocate full fields (different versions)
% Volume, Alternating trapezoidal, Kapur-Rokhlin, and Surface
[B_rV, B_zV, B_rLAlt, B_zLAlt, B_rLKR, B_zLKR, B_rS, B_zS] ...
    = deal( zeros(size(t0t0)) );

for j = 1:2*Mx
    
    fprintf('Full field iter %d/%d.\n', j, 2*Mx)
    
    t0 = t0t0(j);
    
    % compute full fields
    [B_rV(j), B_zV(j)]          =        volume_integral_Jungpyo(RR(j), ZZ(j), shape, 'gauss', 2*My);
    [B_rLAlt(j), B_zLAlt(j)]    = Zakharov_line_integral_Jungpyo(RR(j), ZZ(j), shape, My, true, t0, true);  % line integral on surface at t0; alternating trapezoidal rule
    [B_rS(j), B_zS(j)]          =       surface_integral_Jungpyo(RR(j), ZZ(j), shape, My, true, t0, false); % surface integral on surface at t0
    
    if My >= order+1 % KR quadrature corrections are defined
        [B_rLKR(j), B_zLKR(j)]  = Zakharov_line_integral_Jungpyo(RR(j), ZZ(j), shape, My, true, t0, false, order); % line integral on surface at t0; KR
    end
    
end

% dot with normal

normalnorm = sqrt(dR_dt0.^2 + dZ_dt0.^2);
B_normal_V      = (dZ_dt0.*B_rV     - dR_dt0.*B_zV)     ./ normalnorm;
B_normal_LAlt   = (dZ_dt0.*B_rLAlt  - dR_dt0.*B_zLAlt)  ./ normalnorm;
B_normal_LKR    = (dZ_dt0.*B_rLKR   - dR_dt0.*B_zLKR)   ./ normalnorm;
B_normal_S      = (dZ_dt0.*B_rS     - dR_dt0.*B_zS)     ./ normalnorm;

%% plot different methods of dot(n, B)

figure(1); clf;
plot(t0t0, B_normal_vecpot, 'b-')
hold on
plot(t0t0, B_normal_LKR, 'r--')
plot(t0t0, B_normal_LAlt, 'g-')
plot(t0t0, B_normal_S, 'm-.')
plot(t0t0, B_normal_V, 'cs')
line([t0t0(1), t0t0(end)], [0, 0], 'color', 'black')
xlabel('$t_0$','interpreter','latex')
ylabel('$\vec n \cdot \vec B$','interpreter','latex')
legend('vecpot(KR)', 'line(KR)', 'line(Alt)', 'surface', 'volume')
axis tight
% saveas(gcf, './img/B_normal_is_smooth.png')