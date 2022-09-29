function I = clenshaw_curtis(f,n,a0,b0) % (n+1)-pt C-C quadrature of f
% by L. Trefethen in "Is Gauss Quadrature Better
% than Clenshaw–Curtis?"

x = cos(pi*(0:n)'/n);                       % Chebyshev points
x = ((b0-a0)/2) .* x + (a0+b0)/2;           % conversion to [a,b]
fx = feval(f,x)/(2*n);                      % f evaluated at these points
g = real(fft(fx([1:n+1 n:-1:2])));          % Fast Fourier Transform
a = [g(1); g(2:n)+g(2*n:-1:n+2); g(n+1)];   % Chebyshev coefficients
w = 0*a'; w(1:2:end) = 2./(1-(0:2:n).^2);   % weight vector
I = (b0-a0)/2 * (w*a);                      % the integral