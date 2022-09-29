%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compare my dot(n, B) calculations with those by Dhairya Malhotra
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; clc; close all

%% Parameterize surface

% shape parameters
shape.kappa = 1.7;
shape.a = 1/3;
shape.R0 = 1;
shape.FB = 1;
shape.q0 = 1;

% discretization numbers
Mx = 100;

% observation point discretization
hx = pi/Mx;
t0t0 = hx.*(-Mx : Mx-1)' + pi; % [-pi, pi] + pi

% generate "true" data from a high-resolution self-computation
if isfile('data/B_normal_Kappa1Point7_self_highres.mat')
    load('B_normal_Kappa1Point7_self_highres.mat')
else
    My = 1e3; % number of quadrature nodes per observation point; high resolution
    order = 10;
    B_normal_highres = B_normal_Jungpyo(shape, t0t0, My, order);
    save('data/B_normal_Kappa1Point7_self_highres', 'B_normal_highres')
end

% KR convergence parameter
order = 10;

%% vector potential method

% max My ~ 300 for order 10; 5000 for order 6; 25000 for order 2
MyMy = round( exp(linspace(log(20), log(350), 16)) )';
err  = zeros(size(MyMy));
for j = 1:length(MyMy)
    
    fprintf('Iter %d/%d...', j, length(MyMy))
    My = MyMy(j);
    
    % compute dot(n, B)
    B_normal_vecpot = B_normal_Jungpyo(shape, t0t0, My, order);
    err(j) = norm( B_normal_vecpot - B_normal_highres, inf ); % maximum difference over all observation x's
    
    fprintf('Done.\n')
end
err = err ./ norm(B_normal_highres, inf); % normalize

%% plot comparison with hig-resolution value

figure(1); clf;
loglog(2*MyMy-1, err, 'bx-');
xlabel('Number of quadrature nodes')
ylabel('Normalized maximum difference')
hold on
loglog(2*MyMy(MyMy<200)-1, err(1)*(MyMy(1)./MyMy(MyMy<200)).^(order), 'k:')
legend('Quadrature error', sprintf('$O(h^{%d})$', order),'interpreter','latex')
axis tight
ax = gca;
ax.XAxis.TickValues = sort(unique( [2*MyMy(1)-1, 2*MyMy(end)-1, 100] ));
% saveas(gcf, sprintf('./img/Dhairya_diff_selfconvergence_%dpts_order%d_Kappa1Point7.eps', 2*Mx, order), 'epsc')