function I = gauss(f, n, a, b) % (n+1)-pt Gauss quadrature of f
% by L. Trefethen in "Is Gauss Quadrature Better
% than Clenshaw–Curtis?"

beta = .5./sqrt(1-(2*(1:n)).^(-2));     % 3-term recurrence coeffs
T = diag(beta,1) + diag(beta,-1);       % Jacobi matrix
[V,D] = eig(T);                         % eigenvalue decomposition
x = diag(D); [x,i] = sort(x);           % nodes (= Legendre points) on [-1,1]
w = 2*V(1,i).^2;                        % weights on [-1,1]

x = ((b-a)/2) .* x + (a+b)/2;           % conversion to [a,b]

I = ((b-a)/2)*(w*feval(f,x));           % the integral