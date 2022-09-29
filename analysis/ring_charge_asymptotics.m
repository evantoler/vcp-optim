%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Verify that the potential due to a ring of charges decays at infinity.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% place a ring of charge at (r0, z0) and rotated through toroidal angle
r0 = 2;
z0 = 3;

angle = 1; % pick a direction angle in [-pi/2, pi/2] in the (r,z) plane

RADS = logspace(0, 6, 16); % increasing distances from the origin

magnitudes = nan*RADS; % magnitude of potential at distances in RADS
for j = 1:length(RADS)
    RAD = RADS(j); % RAD = sqrt(r^2 + z^2) = dist to origin
    r = RAD*cos(angle);
    z = RAD*sin(angle);
    
    % evaluate ring charge potential function "u"
    alpha = r^2 + r0^2 + (z-z0)^2;
    beta = 2 * r * r0;
    k2 = 2*beta / (alpha+beta);
    u = 4 * r0 / sqrt(alpha + beta) * ellipticK(k2);

    magnitudes(j) = abs(u);
end

figure;
loglog(RADS, magnitudes, 'bx-')
hold on
loglog(RADS, 1./RADS, 'k:')
legend('Function', '1/Distance')
xlabel 'Distance from origin'
ylabel 'Function magnitude'