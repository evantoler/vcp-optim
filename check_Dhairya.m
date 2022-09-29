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

% read in "true" data from Dhairya's code
do_HighRes = false; % use high-resolution data for kappa=1?
if shape.kappa == 1.7
    B_normal_Dhairya = textread('Dhairya_kappa1point7_Np1200.txt');
    B_normal_Dhairya = B_normal_Dhairya(1:6:end); % get 200 values of t0 instead of 1200
elseif shape.kappa == 1
    if do_HighRes
        B_normal_Dhairya = textread('HighRes_Data.txt');
        B_normal_Dhairya = B_normal_Dhairya(1:6:end); % get 200 values of t0 instead of 1200 for HighRes
    else
        B_normal_Dhairya = textread('Comparison_DataPoints.txt');
    end
else
    error('No comparison data.')
end
B_normal_Dhairya = -B_normal_Dhairya; % align sign conventions
N_Dhairya = length(B_normal_Dhairya);

% discretization numbers
Mx = N_Dhairya/2;

% observation point discretization
hx = pi/Mx;
t0t0 = hx.*(-Mx : Mx-1)' + pi; % [-pi, pi] + pi

% KR convergence parameter
order = 10;

%% vector potential method

% max My ~ 300 for order 10; 5000 for order 6; 25000 for order 2
MyMy = round( exp(linspace(log(20), log(350), 16)) )';
err  = zeros(size(MyMy));

% loop over different numbers of quadrature nodes
for j = length(MyMy):-1:1
    
    fprintf('Iter %d/%d...', j, length(MyMy))
    My = MyMy(j);
    
    % compute dot(n, B) by my code and track discrepancies
    B_normal_vecpot = B_normal_Jungpyo(shape, t0t0, My, order);
    err(j) = norm( B_normal_vecpot - B_normal_Dhairya, inf );
    
    fprintf('Done.\n')
end
err = err ./ norm(B_normal_Dhairya, inf); % normalize

%% Plot comparison with Dhairya values

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
% saveas(gcf, sprintf('./img/Dhairya_diff_convergence_%dpts_order%d.eps', 2*Mx, order), 'epsc')