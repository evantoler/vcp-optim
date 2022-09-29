function [w, beta] = KR_weights(order, varargin)
%{
Return the correction weights from Table 6 of Kapur-Rokhlin.
Corrections are for a logarithmic singularity.

INPUTS:
    order (in {2,6,10}) convergence order; analytically the quadrature
        error decreases as O(h^order), where h = grid spacing.
        Called "k" in the Kapur-Rokhlin paper.
    smooth (in {3,5,9,17,33,43}) interval endpoint correction parameter.
        Must be large enough to match the smoothness of the integrand.
        Called "m" in the Kapur-Rohklin paper.

OUTPUTS:
    w (2*order) correction weights at grid points near the singularity.
        sorted as [-order, ..., -1, 1, ..., order].
        Called "gamma" in the Kapur-Rohklin paper.
    beta((smooth-1)/2) correction weights at new nodes beyond the 
        integration endpoints. Sorted as [1, ..., (smooth-1)/2].
%}


switch order
    case 0
        w = [];
    case 2
        w = [-0.6032109664493744D+00; 0.7518812338640025D+00;
            0.1073866830872157D+01; -0.7225370982867850D+00];
    case 6
        w = [-0.8837770983721025; 0.4799117710681772e1; -0.1064623987147282e2; 0.1219590847580216e2; -0.7407035584542865e1; 0.2051970990601252e1; ...
            0.2915391987686506e1; -0.8797979464048396e1; 0.1365562914252423e2; -0.1157975479644601e2; 0.5130987287355766e1; -0.9342187797694916];
    case 10
        w = [-0.1066655310499552D+01; 0.1009036069527147D+02; -0.4269031893958787D+02; 0.1061953812152787D+03; -0.1715855846429547D+03; 0.1874446431742073D+03; -0.1393153744796911D+03; 0.6872858265408605D+02; -0.2096116396850468D+02; 0.3256353919777872D+01;
            0.4576078100790908D+01; -0.2469045273524281D+02; 0.7648830198138171D+02; -0.1508194558089468D+03; 0.1996415730837827D+03; -0.1807965537141134D+03; 0.1110467735366555D+03; -0.4438764193424203D+02; 0.1044548196545488D+02; -0.1100328792904271D+01];
    otherwise
        error('Input ''order'' must be one of {2,6,10}.')
end

if nargout > 1
    smooth = varargin{1};
    switch smooth
        case 0
            beta = [];
        case 3
            beta = 0.4166666666666667E-01;
        case 5
            beta = [0.5694444444444444E-01; -.7638888888888889E-02];
        case 9
            beta = [0.6965636022927689E-01;
                -.1877177028218695E-01;
                0.3643353174603175E-02;
                -.3440531305114639E-03];
        case 17
            beta = [0.7836226334784645E-01;
                -.2965891540255508E-01;
                0.1100166460634853E-01;
                -.3464763345380610E-02;
                0.8560837610996298E-03;
                -.1531936403942661E-03;
                0.1753039202853559E-04;
                -.9595026156320693E-06];
        case 33
            beta = [0.8356586223906431E-01;
                -.3772568901686391E-01;
                0.1891730418359046E-01;
                -.9296840793733075E-02;
                0.4266725355474016E-02;
                -.1781711570625946E-02;
                0.6648868875120770E-03;
                -.2183589125884841E-03;
                0.6214890604453148E-04;
                -.1506576957395117E-04;
                0.3044582263327824E-05;
                -.4984930776384444E-06;
                0.6348092751221161E-07;
                -.5895566482845523E-08;
                0.3550460453274996E-09;
                -.1040273372883201E-10];
        case 43
            beta = [0.8490582345073516E-01;
                -.4001723785254229E-01;
                0.2156339227395192E-01;
                -.1173947578371037E-01;
                0.6165108551649839E-02;
                -.3051271143145484E-02;
                0.1403005122150106E-02;
                -.5931791433462842E-03;
                0.2286250628123645E-03;
                -.7968542809070158E-04;
                0.2490991825767152E-04;
                -.6921164516465828E-05;
                0.1691476513287747E-05;
                -.3590633248885163E-06;
                0.6517156577922871E-07;
                -.9908863655077215E-08;
                0.1227209060809220E-08;
                -.1188834746888414E-09;
                0.8447408532519018E-11;
                -.3914655644778233E-12;
                0.8806394737861057E-14];
        otherwise
            error('Input value for ''smooth'' must be one of {3,5,9,17,33,43}.')
    end
end