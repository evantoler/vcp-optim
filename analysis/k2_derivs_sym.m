clc; clear;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Determine the behavior of a coefficient that appears in a numerical test
% for harmonic functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

syms r(t) z(t) R Z t0

% % outdated coefficient test
% k2 = 4 * R * r(t) / ( (R + r(t))^2 + (Z - z(t))^2 );
% dk2dR = diff(k2, R);
% dk2dZ = diff(k2, Z);

% pretty(simplify(dk2dR))
% pretty(simplify(dk2dZ))

% coef1 = (R+r)*(1-2/k2) + 4*r*(1 - 2*R*(R+r)/((R + r(t))^2 + (Z - z(t))^2)) / k2 * (-2/k2 + 1/2);
% % pretty(simplify(coef1))
% coef2 = (R+r)*2/k2 + 4*r*(1 - 2*R*(R+r)/((R + r(t))^2 + (Z - z(t))^2)) / k2 ...
%     * ( 2/k2 + 1/k2*(1/(1-k2)-1) - 1/(2*k2) );
% pretty(simplify(coef2))

% k2 = 4 * r(t0) * r(t) / ( (r(t0) + r(t))^2 + (z(t0) - z(t))^2 );

%% coefficient in harmonic test
coef = diff(z(t)) * (r(t0) - r(t)) - diff(r(t)) * (z(t0) - z(t));
coef = coef / ( (r(t0) - r(t))^2 + (z(t0) - z(t))^2 );

%% check derivatives for L'Hopital rule
coef1 = simplify( diff(coef, t) );
% L'hopital - numerator and denomenator
coef1num = (((r(t) - r(t0))*diff(z(t), t) - (z(t) - z(t0))*diff(r(t), t))*(2*(r(t) - r(t0))*diff(r(t), t) + 2*(z(t) - z(t0))*diff(z(t), t))) - ((r(t) - r(t0))*diff(z(t), t, t) - (z(t) - z(t0))*diff(r(t), t, t))*((r(t) - r(t0))^2 + (z(t) - z(t0))^2);
coef1den = ((r(t) - r(t0))^2 + (z(t) - z(t0))^2)^2;
dcoef1num = diff(coef1num, t, t, t, t);     % view iterated derivatives
dcoef1den = diff(coef1den, t, t, t, t);
limitdcoef = dcoef1num / dcoef1den;
pretty( simplify(subs(limitdcoef, t, t0)) )