clear; clc;

% definition of modulus k^2 for points on the generating curve
syms r(t) z(t) t0
k2 = 4 * r(t0) * r(t) / ( (r(t0) + r(t))^2 + (z(t0) - z(t))^2 );

% take analytical derivatives
dk2 = simplify( diff(k2, t) );      % first derivative
d2k2 = simplify( diff(dk2, t) );    % second derivative
d3k2 = simplify( diff(d2k2, t) );   % third derivative

disp('first deriv')
dk2 = simplify( subs(dk2, t, t0) ); pretty(dk2)
disp('second deriv')
d2k2 = simplify( subs(d2k2, t, t0) ); pretty(d2k2)
disp('third deriv')
d3k2 = simplify( subs(d3k2, t, t0) ); pretty(d3k2)