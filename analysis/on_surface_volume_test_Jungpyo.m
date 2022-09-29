%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% test asymptotic behavior of a principal value quantity
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; close all; clc

% set shape parameters
kappa = 1;
a = 1/3;
R0 = 1;
FB = 1;
q0 = 1;

Rmax = sqrt( R0^2 + 2*a*R0 );
zmax = @(r) kappa*a*R0./r .* sqrt(1 - ((r.^2-R0^2)./(2*a*R0)).^2);

R = Rmax;
Z = 0;

r = R;

% elliptic integral parameters
k2 = 4*r*R ./ ( (r+R)^2 + (zmax(r))^2 );
n = 4*r*R ./ (r+R)^2;

% quantity to test
f = @(epsi) ellipticF((pi-epsi)./2,k2) + ...
    (r-R)./(r+R) .* ellipticPi(n, (pi-epsi)./2, k2);

figure;
low = 1e-3;
fplot(f, [low, 1]); hold on
fplot(@(x) (low*f(low))./x, [low, 1])
fplot(@(x) -log(x), [low, 1])
xlabel('x')
legend('function','O(1/x)','O(|log(x)|)')