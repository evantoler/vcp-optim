function [B_r, B_z] = volume_integral_Jungpyo(R, Z, shape, varargin)
%{
INPUTS:
R,Z target point
shape = structure of parameters
N+1 quadrature points
varargin:
    method: Integration method; options are 'gauss', 'kr', 'builtin'
    N       number of quadrature nodes for 'gauss' and 'kr' quadratures

OUTPUTS:
    B_r     (scalar) radial component of the magnetic field at (R,Z)
    B_z     (scalar) vertical """"""""""""""""""""""""""""""""""""""

Evan Toler, 2022
%}

if nargin>3, method=varargin{1}; else, method='builtin'; end
method = lower(method); % conform to lowercase
if ~strcmpi(method,'builtin') 
    if nargin>4, N=varargin{2}; else, N=100; end 
end

% surface shape parameters
kappa = shape.kappa;
a  = shape.a;
R0 = shape.R0;
FB = shape.FB;
q0 = shape.q0;
Cs = FB*(kappa+1/kappa)/(R0^3*q0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% bounds on coordinates
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Rmin = sqrt( R0^2 - 2*a*R0 );
Rmax = sqrt( R0^2 + 2*a*R0 );
zmax = @(r) kappa*a*R0./r .* sqrt(1 - ((r.^2-R0^2)./(2*a*R0)).^2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% my reduction to a 1D integral
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% important quantities
k2_plus  = @(r) 4*R.*r ./ ( (R + r).^2 + (Z + zmax(r)).^2 );
k2_minus = @(r) 4*R.*r ./ ( (R + r).^2 + (Z - zmax(r)).^2 );
n = @(r) 4*R.*r ./ (R+r).^2;

% build integrands in each coordinate
plus_term  = @(r) 1./sqrt(k2_plus(r)) ...
    .* ( (1-k2_plus(r)./2).*ellipticK(k2_plus(r))...
    - ellipticE(k2_plus(r)) );
minus_term = @(r) 1./sqrt(k2_minus(r)) ...
    .* ( (1-k2_minus(r)./2).*ellipticK(k2_minus(r))...
    - ellipticE(k2_minus(r)) );
B_r = @(r) Cs/pi .* r.^(3/2) ./ sqrt(R) ...
    .* (plus_term(r) - minus_term(r));

plus_term  = @(r) (Z + zmax(r)) ./ sqrt( (R+r).^2 + (Z+zmax(r)).^2 ) ...
    .* ( ellipticK(k2_plus(r)) ...
    + (r-R)./(r+R) .* ellipticPi(n(r), k2_plus(r)) );
minus_term = @(r) (Z - zmax(r)) ./ sqrt( (R+r).^2 + (Z-zmax(r)).^2 ) ...
    .* ( ellipticK(k2_minus(r)) ...
    + (r-R)./(r+R) .* ellipticPi(n(r), k2_minus(r)) );
B_z = @(r) real( ...
    -Cs/(2*pi) .* r ...
    .* (plus_term(r) - minus_term(r)) ...
    );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot integrand
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% figure; 
% subplot(1,2,1)
% fplot(B_r, [Rmin, Rmax]) % plot integrand
% title('Integrand, $B_r = \int_{R_{min}}^{R_{max}} f(r) \, dr$', 'interpreter', 'latex')
% xlabel $r$
% ylabel $f(r)$
% 
% subplot(1,2,2)
% fplot(B_z, [Rmin, Rmax]) % plot integrand
% title('Integrand, $B_z = \int_{R_{min}}^{R_{max}} f(r) \, dr$', 'interpreter', 'latex')
% xlabel $r$
% ylabel $f(r)$

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% integrate
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch method
    
    case 'gauss' % Gauss quadrature
        
        B_r = gauss(B_r, N, Rmin, Rmax);
        B_z = gauss(B_z, N, Rmin, Rmax);
        
    case 'kr' % Kapur-Rokhlin quadrature (aperiodic version)
        
        % get KR parameters
        order = 6;
        smooth = 3;
        [gamma, beta] = KR_weights(order, smooth); 
        
        % equispaced grid in radius
        rr = linspace(Rmin, Rmax, N);
        h = rr(2)-rr(1);
        idx=( 1:(smooth-1)/2 )';
        
        % compute integral w/ corrections
        B_r = h*B_r(rr(1))/2 + h*sum( B_r(rr(2:N-1)) ) ...
            + h*dot(beta, B_r(rr(1 + idx)) - B_r(Rmin - h*idx) ) ... % left (smooth) correction
            + h*dot(gamma, B_r(Rmax+h*[order:-1:1, -1:-1:-order])); % right (singular) correction
        
        B_z = h*B_z(rr(1))/2 + h*sum( B_z(rr(2:N-1)) ) ...
            + h*dot(beta, B_z(rr(1 + idx)) - B_z(Rmin - h*idx) ) ... % left (smooth) correction
            + h*dot(gamma, B_z(Rmax+h*[order:-1:1, -1:-1:-order])); % right (singular) correction
        
    otherwise % builtin
        
        B_r = integral(B_r, Rmin, Rmax, 'abstol', 1e-15);
        B_z = integral(B_z, Rmin, Rmax, 'abstol', 1e-15);
        
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% outdated triple integral version
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% denom = @(r,theta,z) ( R^2 + r.^2 + (Z-z).^2 - 2*R*r.*cos(theta) ).^(3/2);
% B_r   = @(r,theta,z) r.^2 .* (Z-z) ./ denom(r,theta,z);
% B_z   = @(r,theta,z) r.^2 .* (r - R*cos(theta)) ./ denom(r,theta,z);
% 
% % zmin = @(r) -zmax(r);
% B_r = integral3(B_r, Rmin, Rmax, -pi, 0, @(r,theta)-zmax(r), @(r,theta)zmax(r), 'abstol', 1e-10, 'method', 'iterated') ...
%     + integral3(B_r, Rmin, Rmax, 0, pi,  @(r,theta)-zmax(r), @(r,theta)zmax(r), 'abstol', 1e-10, 'method', 'iterated');
% B_z = integral3(B_z, Rmin, Rmax, -pi, 0, @(r,theta)-zmax(r), @(r,theta)zmax(r), 'abstol', 1e-10, 'method', 'iterated') ...
%     + integral3(B_z, Rmin, Rmax, 0, pi,  @(r,theta)-zmax(r), @(r,theta)zmax(r), 'abstol', 1e-10, 'method', 'iterated');

return