function [ elem ] = OrbitElem( R, V, mu )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% OrbitElem 
%  |Converts cartesian orbital elements to classical|
%   Ouput: [p,a,e,i,RAAN,w,v,E,M]
%           Special cases:
%               v - true longitude of periapsis for non-circular equatorial
%            orbits
%               v - argument of latitude for circular inclined orbit
%               v - True longitude for circular Equatorial orbit   
%   Input: R=[x,y,z], V=[xdot, ydot, zdot]
%   Units: rad, km, sec 
%   Tolerance: 1e-17
%   Reference: Vallado
% Author: Bharat Mahajan (https://github.com/princemahajan)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% CONSTANTS
TOL = 1e-15;

r = norm(R);
v = norm(V);
rinv = 1/r;

H = cross(R,V); % specific momentum in cartesian system
h = norm(H);

n_vec = cross([0;0;1],H); % Line of nodes vector in cartesian
n = norm(n_vec);
ninv = 1/n;
if n < TOL
    % Zero inclination, so node line coincides with X-axis
    n_vec = [1;0;0];
    ninv = 1/norm(n_vec);
end
    
% compute eccentricity vector
e_vec = cross(V,H)/mu - R*rinv;
e = norm(e_vec);
einv = 1/e;
if e < TOL
    % Zero eccentricity, so eccentricity vector coincides with node line
    e_vec = n_vec;
    einv = 1/norm(e_vec);
end

p = (h^2)/mu; % semi-latus rectum

% SMA is not defined for parabolic orbits
if abs(1-e) < TOL
    a = Inf;
else
    a = -mu/2/(v^2/2 - mu/r); % Vis-viva equation
end

i = acos(H(3)/h); % inclination between 0 and pi

raan = acos(n_vec(1)*ninv); % RAAN
% do quad check
if (n_vec(2) < 0)
    raan = 2*pi-raan;
end

w = acos(dot(n_vec,e_vec)*ninv*einv); % Argument of Periapsis
% do quad check
if (e_vec(3) < 0)
    w = 2*pi-w;
end
if abs(imag(w)) < 1e-7
    w = real(w);
else
    disp('OrbitElem: imaginary AOP');
    return;
end

% True anomaly
cv = dot(e_vec,R)*einv*rinv;
if abs(abs(cv) - 1) < TOL
    cv = sign(cv);
end

v = acos(cv); % True Anamoly
% do quad check
rvp = dot(R,V);
if (rvp < 0)
    v = 2*pi-v;
end


% Special Cases
if (e >= TOL) && (i < TOL || abs(i-pi) < TOL)
    % Elliptical Equatorial
    if (i == pi) && (e_vec(2) > 0)
        % retrograde motion
        w = 2*pi-w;
    elseif (e_vec(2) < 0)
        w = 2*pi-w;
    end
elseif (i >= TOL) && (e < TOL)
    % Circularly Inclined (no periapsis)
    v = acos(dot(n_vec,R)*rinv*ninv);
    % do quad check
    if (R(3) < 0)
        v = 2*pi-v;
    end
elseif (i < TOL) && (e < TOL)
    % Circular Equatorial orbit
    v = acos(R(1)*rinv);
    % do quad check
    if (R(2) < 0)
        v = 2*pi-v;
    end
end

if a > 0
    % find M and E
    E = 2*atan(sqrt((1-e)/(1+e))*tan(v/2));
    if E < 0
        E = 2*pi + E;
    end
    M = E - (e*sin(E));
    if M < 0
        M = 2*pi + M;
    end
else
    E = 0;
    M = 0;
end

% Prepare Output
if abs(i) < TOL
    raan = 0;
end
if abs(e) < TOL
    w = 0;
end

elem = [p,a,e,i,raan,w,v,E,M];

end

