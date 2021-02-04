function Xdot = GravAcc(t,X_GCRF,t0,GMST0, we, degree, order, mu, Re, Clm, Slm)
%GravAcc Computes acceleration due to Zonal and Tesseral Gravity Harmonics
%in the ECI frame. Max 70x70.
%
% Clm, Slm: Array of C, S unnormalized coefficients
% we: Earth rotational speed in radians/s
%
% Author: Bharat Mahajan (https://github.com/princemahajan)


MAX_DEGREE = 70;
MAX_ORDER = 70;

% Make sure order is less than degree
if degree <= MAX_DEGREE && order <= MAX_ORDER && order <= degree
    l_max = degree;
    m_max = order;
elseif degree > MAX_DEGREE && order <= MAX_ORDER && order <= MAX_DEGREE
    l_max = MAX_DEGREE;
    m_max = order;
    disp(['Degree or order changed to ' num2str(l_max) ',' num2str(m_max)]);
elseif degree <= MAX_DEGREE && (order >= MAX_ORDER || order >= degree)
    l_max = degree;
    m_max = degree;
    disp(['Degree or order changed to ' num2str(l_max) ',' num2str(m_max)]);
else
    l_max = MAX_DEGREE;
    m_max = MAX_ORDER;
    disp(['Degree or order changed to ' num2str(l_max) ',' num2str(m_max)]);
end

% disp(['Gravity Spherical Expansions Degree, Order: ' num2str(l_max) ', ' num2str(m_max)]);

% Current Greenwhic mean sidereal time or Earth Rotation Angle
GMST = GMST0 + we*(t - t0);

% DCM between ECI and ECF frames
R3 = @(th) [cos(th), sin(th), 0;
            -sin(th), cos(th), 0;
            0,        0,       1];

DCM = R3(GMST); 

% Position vector in Terrestial frame ECF
X_ITRF = DCM*X_GCRF(1:3,1);
R = X_ITRF';

% radial distance
r = norm(X_ITRF);

% latitude between -pi/2 and pi/2
lat = asin(X_ITRF(3)/r);

% longitude between 0 and 2*pi
lon = atan2(X_ITRF(2),X_ITRF(1));

% Compute associated legendre polynomials
PlmM = zeros(l_max+1, l_max+2);
for l = 2:l_max
    PlmM(l,1:l+1) = legendre(l,sin(lat))';
end

% correct for the negative sign in the Plm definition of Matlab
Plm = @(l,m) PlmM(l,m)/(-1)^(m-1);

% Compute nonspherical potential partials (See Vallado page 550)
dUdr = 0;
dUdlat = 0;
dUdlon = 0;

for l = 2:l_max

    for m = 0:min(m_max,l)
   
        CS1 = Clm(l+1,m+1)*cos(m*lon) + Slm(l+1,m+1)*sin(m*lon);
        CS2 = Slm(l+1,m+1)*cos(m*lon) - Clm(l+1,m+1)*sin(m*lon);
     
        % dUdr
        dUdrlm = -mu/r^2*(Re/r)^l*(l+1)*Plm(l,m+1)*CS1;

        % dUdlat
        dUdlatlm = mu/r*(Re/r)^l*(Plm(l,m+2) - m*tan(lat)*Plm(l,m+1))*CS1;
      
        % dUdlon
        dUdlonlm = mu/r*(Re/r)^l*m*Plm(l,m+1)*CS2;

        dUdr = dUdr + dUdrlm;
        dUdlat = dUdlat + dUdlatlm;
        dUdlon = dUdlon + dUdlonlm; 
    end
end


% Acceleration in ITRF frame
rI = R(1);
rJ = R(2);
rK = R(3);
rIJ = sqrt(R(1)^2 + R(2)^2);

asph_ITRF = -mu/r^3*R';

aI_ITRF = ((1/r*dUdr - rK/(r^2*rIJ)*dUdlat)*rI - 1/rIJ^2*dUdlon*rJ);

aJ_ITRF = ((1/r*dUdr - rK/(r^2*rIJ)*dUdlat)*rJ + 1/rIJ^2*dUdlon*rI);
            
aK_ITRF = (1/r*dUdr*rK + rIJ/r^2*dUdlat);
  
% aIJK_ITRF = 1.5*Clm(3,1)*(mu/r^2)*(Re/r)^2*[(1-5*(rK/r)^2)*(rI/r); (1-5*(rK/r)^2)*(rJ/r); (3-5*(rK/r)^2)*(rK/r)];

% disp(aI_ITRF - aIJK_ITRF(1));
% disp(aJ_ITRF - aIJK_ITRF(2));

a_ITRF = asph_ITRF + [aI_ITRF; aJ_ITRF; aK_ITRF];

% Transform to GCRF
Xdot(1:3,1) = X_GCRF(4:6,1);
Xdot(4:6,1) = DCM'*a_ITRF;

end



