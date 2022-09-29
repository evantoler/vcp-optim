close all; clear; clc

shape.kappa = 1;
shape.a = 1/3;
shape.R0 = 1;
shape.FB = 1;
shape.q0 = 1;

% an observation point outside the generating curve
R = 1.5;
Z = 0.4;

% generate high-precision "true" field values
if exist('Bstar_Jungpyo.mat', 'file')
    load Bstar_Jungpyo;
else
    N = 20000;
    
    [B_rV, B_zV] = volume_integral_Jungpyo(R, Z, shape, N);
    [B_rZ, B_zZ] = Zakharov_line_integral_Jungpyo(R, Z, shape, N);
    
    % verify that both computations are valid
    disp 'Relative B field difference:'
    disp(norm([B_rZ; B_zZ] - [B_rV; B_zV])/norm([B_rZ; B_zZ]))
    
    B_rSTAR = B_rV; B_zSTAR = B_zV;
    save('Bstar_Jungpyo.mat', 'B_rSTAR', 'B_zSTAR');
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% convergence and time scaling with number of quadrature points
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

norm_Bstar = sqrt(B_rSTAR^2 + B_zSTAR^2);

% NN = round(exp(linspace(log(5), log(200), 20)))';
NN = round(linspace(5, 100, 20));

% preallocate
[errV, timeV, errZ, timeZ] ...
    = deal( zeros(size(NN)) ); % V=volume, Z=Zakharov

for j = length(NN):-1:1
    
    fprintf('Iter%d\n',j)
    N = NN(j);
    
    % volume integral calculation
    tic;
    [B_rV, B_zV] = volume_integral_Jungpyo(R, Z, shape, 'gauss', N);
    timeV(j) = toc;
    errV(j) = norm([B_rV; B_zV] - [B_rSTAR; B_zSTAR])/norm_Bstar;
    
    % line integral calculation
    tic;
    [B_rZ, B_zZ] = Zakharov_line_integral_Jungpyo(R, Z, shape, N);
    timeZ(j) = toc;
    errZ(j) = norm([B_rZ; B_zZ] - [B_rSTAR; B_zSTAR])/norm_Bstar;
    
end

%% Plotting

figure;
semilogy(NN, errV, 'bx-', NN, errZ, 'ro-'); hold on
semilogy(NN, 10*NN.^(-3), 'k--')
semilogy(NN(NN<40), exp(-NN(NN<40)), 'k:')
xlabel('$N$ (\# quadrature points)','fontsize',20,'interpreter','latex')
ylabel('$\|\vec B_{approx} - \vec B_*\|_2/\|\vec B_*\|_2$','fontsize',20,'interpreter','latex')
title('Jungpyo''s geometry','fontsize',20)
axis tight
legend('volume', 'Zakharov','$O(N^{-3})$', '$O(e^{-N})$', ...
    'interpreter','latex', 'location', 'east')
% saveas(gcf, 'error_convergence_Zakharov_as_trapezoid_Jungpyo.jpg')

figure;
semilogy(NN, timeV, 'bx-', NN, timeZ, 'ro-')
xlabel('$N$ (\# quadrature points)','fontsize',20,'interpreter','latex')
ylabel('time (sec)','fontsize',20)
title('Jungpyo''s geometry','fontsize',20)
axis tight
legend('volume', 'Zakharov', 'location', 'northwest')

