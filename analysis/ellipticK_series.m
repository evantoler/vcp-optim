function val = ellipticK_series(k2)
%{
Evaluate the complete elliptic integral
    \int_0^{\pi/2} d \theta / \sqrt(1 - k^2 \sin^2 \theta)
by a series expansion valid when the modulus k2 ("k squared") 
is close to 1.

Reference: Gradshteyn & Ryzhik, 8.113.2

Evan Toler, 2022
%}

tol = 1e-13;

kp = sqrt(1-k2);                % complementary modulus "k prime"
m2 = ( (1-kp) ./ (1+kp) ).^2;   % new parameter for series

% sum the series
val = ones(size(k2)); term = val;
j = 1;
while norm(term, inf) > tol
    term = (j/(j+1))^2 .* m2 .* term;
    val = val + term;
    j = j + 2;
end

% leading factor
val = pi./(1+kp) .* val;