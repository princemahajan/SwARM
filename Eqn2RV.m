%% Equinoctial to Cartesian Elements Converter
%   Eqn: [a,ML,p1(=tan(i/2)cos(h)),p2(=tan(i/2)sin(h)),q1(=e cos(g+h)),q2(=e sin(g+h))];
% Author: Bharat Mahajan (https://github.com/princemahajan)

function out = Eqn2RV(in ,mu, tol, reverse)

if reverse == true
   
    R = in(1:3);
    V = in(4:6);
    
    % Classical elements
    [ elem ] = OrbitElem( R, V, mu );
    oe = elem(2:7);
    oe(end) = elem(end);
   
    % if NaN or inf then replace with 0
    for ctr = 1:6
        if isnan(oe(ctr)) == true || isinf(oe(ctr)) == true
            oe(ctr) = 0;
        end
    end
    
    a = oe(1);
    ML = mod(oe(6) + oe(5) + oe(4),2*pi);
    
    p1 = tan(oe(3)/2)*cos(oe(4));
    p2 = tan(oe(3)/2)*sin(oe(4));
    
    q1 = oe(2)*cos(oe(4) + oe(5));
    q2 = oe(2)*sin(oe(4) + oe(5));
    
    out = [a,ML,p1,p2,q1,q2]';
    
else
    
    % convert to Position Velocity in ECI
    a = in(1); Lambda = in(2); p1 = in(3); p2 = in(4); q1 = in(5); q2 = in(6);

    % True Longitude
    [TL, EL] = Mean2TrueLong(Lambda, q1,q2, tol);
    
    % radial distance
    r = a*(1 - q1*cos(EL) - q2*sin(EL));
    
    p = a*(1 - q1^2 - q2^2);
    
    % R,V in Equinoctial frame
    Req = r*[cos(TL); sin(TL); 0];
    Veq = sqrt(mu/p)*[-q2 - sin(TL); q1 + cos(TL); 0];
    
    % DCM from equinoctial to ECI
    DCM = 1/(1+p1^2+p2^2)*[1-p2^2+p1^2, 2*p1*p2, 2*p2;
                         2*p1*p2,   1+p2^2-p1^2, -2*p1;
                         -2*p2,     2*p1,     1-p1^2-p2^2];
                     
   R = DCM*Req;
   V = DCM*Veq;
   
   out = [R', V']';
   
end

end
