%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% integral from -pi to pi of cos^2(t) * log(|sin(t)|/2) dt
% which is 2pi-periodic and singular at the origin
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

I = -pi/4 * (1 + log(16)); % analytic solution
order = 10;

%% Convergence test for quadrature.
% Note: we need at least "order+2" points 
% for the quadrature rule to be defined.
% We use an even number of quadrature nodes (omitting 0) for convenience.
NN = 2*round( exp(linspace(log((order+2)), log(1000/2), 16)) )';
err = zeros(size(NN));
for j = 1:length(NN)
    N = NN(j);
    h = 2*pi/N;
    tt = h * (-N/2+1 : N/2)';
    ff = cos(tt).^2 .* log(abs(sin(tt./2)));
    
    % apply Kapur-Rokhlin corrected quadrature
    IKR = periodic_KR(h, ff, order);
    err(j) = abs(I - IKR);
end

%% plot convergence
figure(2); clf;
loglog(NN-1, err, 'bx-'); hold on
loglog(NN-1, err(1)*(NN(1)./NN).^order, 'k:')
xlabel('Number of quadrature nodes')
ylabel('Quadrature error')
legend('Error', sprintf('Order %d', order))
