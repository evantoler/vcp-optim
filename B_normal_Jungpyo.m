function [B_normal] =  B_normal_Jungpyo(shape, t0t0, My, varargin)
%{
Compute dot(n, B) at points on the surface corresponding to t0t0.
Each inner quadrature uses 2*My points.
Uses the vector potential as an intermediary.

INPUTS
    shape   (struct) Jungypyo's geometry
    t0t0    (2*Mx)   observation point parameters on the generating curve
    My      (scalar) 2*My points for each inner quadrature 
    order   (in {2,6,10}) Kapur-Rokhlin convergence order
    
OUTPUT
    B_normal    (2*Mx) vector of normal components of B
                    B_normal(j) corresponds to t0 = t0t0(j)

Evan Toler, 2022
%}

C = 1/(4*pi);

% KR parameter
if nargin < 4, order = 10; else, order = varargin{1}; end

% outer loop: over x on the surface; loop over t0 in t0t0
[RR, ~, ~, ~, dR_dt0, dZ_dt0] = get_functions_Jungpyo(shape, 0, t0t0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% inner quadrature
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% make sure KR corrections are defined
if My <= order, error('Error. Need My>=%d.', order+1); end

% inner discretization for the generating curve
hy = pi/My;

% preallocate
% toroidal vector potential
A_phi_KR = zeros( size(t0t0) );

% loop over observation points on surface
parfor j = 1:length(t0t0)
    
    t0 = t0t0(j);
    
    tt = hy.*(-My+1 : My)' + t0; %(-pi+h), ..., pi shifted by t0
    
    % set up integrand
    [~, ~, ~, ~, ~, ~, dpsi_dn_bracket, alpha, beta, k2] ...
        = get_functions_Jungpyo(shape, t0, tt);
    K = ellipticK(k2);
    E = ellipticE(k2);
    A_integrand = C * dpsi_dn_bracket .* 4 ./ sqrt(alpha + beta) ...
        .* ( 2./k2.*( K - E ) - K );

    % Kapur Rokhlin quadrature
    A_phi_KR(j) = periodic_KR(hy, A_integrand, order);
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% take tangential derivative
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Mx = length(t0t0)/2;

psi_VC = RR .* A_phi_KR; % poloidal flux function

% Fourier differentiation in t0
B_normal = ifft(1i*[0:Mx-1, 0, -Mx+1:-1]'.*fft(psi_VC));
B_normal = B_normal ./ sqrt(dR_dt0.^2 + dZ_dt0.^2);     % normalization
B_normal = -B_normal ./ RR;                             % factor -1/R from triple product sign change and grad(e_phi)=(e_phi)/R

return