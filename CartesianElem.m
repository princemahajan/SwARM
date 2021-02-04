%% CartesianElem
%   Converts orbital elements to cartesian elements
%   Input: [p,e,i,RAAN,w,M,nu] (provide M or nu, if both there then M will be used)
%           Special cases:
%               v - true longitude of periapsis for non-circular equatorial
%            orbits
%               v - argument of latitude for circular inclined orbit
%               v - True longitude for circular Equatorial orbit   
%   Output: R=[x,y,z], V=[xdot, ydot, zdot]
%   Units: rad, km, sec 
%   Tolerance: 1e-10
%   Reference: Vallado
% Author: Bharat Mahajan (https://github.com/princemahajan)


function [ R, V ] = CartesianElem( orbelem, mu, isTA )

% Constants
TOL = 1e-10;
TOL_M_E = 1e-10;

% get the orbital elements
p=orbelem(1); e=orbelem(2); i=orbelem(3); raan=orbelem(4);
w=orbelem(5); M=orbelem(6); 

% check whether we have nu or M
if isTA == true
    % we have nu
    v = orbelem(6);
else
%     % we have M now compute E and then nu
%     E = fzero(@(x) x-e*sin(x)-M, M, optimset('TolX',TOL_M_E));
% %     E = fzero(@(x) M-x+e*sin(x), M);
%     % find nu using half angle formula
%     v = 2*atan(sqrt((1+e)/(1-e))*tan(E/2));
%     % adjust v value to be within 0 and 2*pi
%     if v < 0
%         v = v + 2*pi;
%     end
[E, v] = KeplerEqSolver(M, e, TOL);
end

% Handle Special cases
if (e < TOL) && (i < TOL) % Circular Equatorial
    w = 0;
    raan = 0;
elseif (e < TOL) && (i >= TOL) % Circular Inclined
    w = 0;
elseif (i < TOL) && (e >= TOL) % Elliptical Equatorial
    w = raan + w;
    raan = 0;

end

% find Perifocal system coordinates PQW
cv = cos(v);
sv = sin(v);
rdi = 1/(1 + e*cv);
smp = sqrt(mu/p);
r = [ p*cv*rdi;
      p*sv*rdi;
      0];
v = [ -smp*sv;
      smp*(e+cv);
      0];

% compute rotation matrix
cr = cos(raan);
cw = cos(w);
sr = sin(raan);
sw = sin(w);
ci = cos(i);
si = sin(i);
rot=[cr*cw-sr*sw*ci, -cr*sw-sr*cw*ci, sr*si;
     sr*cw+cr*sw*ci, -sr*sw+cr*cw*ci, -cr*si;
     sw*si, cw*si, ci];

% convert spherical coordinates to cartesian
R = rot*r;
V = rot*v;

end

