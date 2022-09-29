%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check pencil-and-paper derivatives of the ring charge potential by 
% finite differences
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; clc; close all;

% pick a random location for the ring of charges: (r0,z0) x [0,2pi]
r0 = 10*rand;
z0 = 10*randn;

% pick evaluation points that aren't on the ring charge
r = r0+3;
z = z0-2;

dir = randn(2,1);   % direction for the finite difference

% "base" function value at the evaluation point
alpha = r^2 + r0^2 + (z-z0)^2;
beta = 2 * r * r0;
k2 = 2*beta / (alpha+beta);
u = 4 * r0 / sqrt(alpha + beta) * ellipticK(k2);

% true analytic gradient components
dudr = 4*r0/(alpha+beta)^(3/2) * ( -(r+r0)*ellipticK(k2) ... 
            + 1/(2*r) * (ellipticE(k2)/(1-k2) - ellipticK(k2)) ...
            * ((z-z0)^2 + r0^2 - r^2) );
dudz = -4*r0*(z-z0)/(alpha+beta)^(3/2) * ellipticE(k2)/(1-k2);

hh = logspace(0, -12, 13); % FD step sizes
err = zeros(size(hh));
for j=1:length(hh)

    h = hh(j);

    % finite difference derivative approximation
    rh = r + h*dir(1);
    zh = z + h*dir(2);

    alphah = rh^2 + r0^2 + (zh-z0)^2;
    betah = 2 * rh * r0;
    k2h = 2*betah / (alphah+betah);
    uh = 4 * r0 / sqrt(alphah + betah) * ellipticK(k2h);

    dirgradh = (uh - u)/h;

    % analytic expression for the "dir"-directional derivative
    dirgrad = dir(1)*dudr + dir(2)*dudz;

    err(j) = abs(dirgrad - dirgradh);
end

figure;
loglog(hh, err, 'bx-')
hold on
loglog(hh, hh, 'k:')
legend('FD derivative', 'first order')