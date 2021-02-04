%% DiffEq2Curv 
%   Converts differential equinoctial elements 
%   (a,Mean Longitude,tan(i/2)*cos(h),tan(i/2)*sin(h),e*cos(g+h),e*sin(g+h)) 
%   to the relative position and velocity in the curvilinear frame.
%   Angular velocity of the Chief is computed using the RAAN rate that
%   includes the effects of all the zonals and tessteals provided.
%
% TBD: Check, velocity transformation might have bugs!
%
% Author: Bharat Mahajan (https://github.com/princemahajan)

function Sigma = DiffEq2Curv(X, theta, degree, order, DiffTA, mu, Re, Clm,Slm, tol) 

% Equinoctial
a = X(1); Lambda = X(2); p1 = X(3); p2 = X(4);q1 = X(5); q2 = X(6);

order = min([degree, order]);

% compute true longitude
DelX = Del2Eqn(X,mu, true);
l = DelX(1);g = DelX(2);h = DelX(3);
e = sqrt(q1^2 + q2^2);
psi = Mean2TrueLong(Lambda, q1, q2, tol);
lop = mod(atan2(q2,q1),2*pi); 
f = mod(psi - lop,2*pi);
eta = sqrt(1 - e^2);

% Intermediate Variables
p = a * (-q1 ^ 2 - q2 ^ 2 + 1);
r = p/(1 + q1*cos(psi) + q2*sin(psi)); 
Vr = sqrt((mu / a / (-q1 ^ 2 - q2 ^ 2 + 1))) * (sin(psi) * q1 - q2 * cos(psi));
Vt = sqrt((mu / a / (-q1 ^ 2 - q2 ^ 2 + 1))) * (cos(psi) * q1 + sin(psi) * q2 + 0.1e1);
sigma1 = p1^2 + p2^2;

% dTLML = -(cos(psi) ^ 2 * eta ^ 2 + 0.2e1 * cos(psi) ^ 2 * q1 ^ 2 + 0.2e1 * sin(psi) * cos(psi) * q1 * q2 - cos(psi) ^ 2 + 0.2e1 * cos(psi) * q1 + 0.2e1 * sin(psi) * q2 - eta ^ 2 - q1 ^ 2 + 0.2e1) / eta ^ 3;
dTLML = 0.1e1 / eta ^ 3 * (0.1e1 + e * cos(f)) ^ 2;
dTLq1 = -(-0.2e1 * cos(psi) ^ 2 * q1 ^ 2 * q2 + cos(psi) * sin(psi) * eta ^ 2 * q1 + 0.2e1 * cos(psi) * sin(psi) * q1 ^ 3 + cos(psi) ^ 2 * eta * q2 - cos(psi) * sin(psi) * eta * q1 + cos(psi) ^ 2 * q2 - 0.2e1 * sin(psi) * cos(psi) * q1 - 0.2e1 * cos(psi) * q1 * q2 + 0.2e1 * sin(psi) * q1 ^ 2 - eta ^ 2 * q2 + q1 ^ 2 * q2 - 0.2e1 * sin(psi) * eta - 0.2e1 * eta * q2 - 0.2e1 * sin(psi) - 0.2e1 * q2) / (0.1e1 + eta) / eta ^ 3;
dTLq2 = -(0.2e1 * cos(psi) ^ 2 * eta ^ 2 * q1 + 0.2e1 * cos(psi) ^ 2 * q1 ^ 3 + cos(psi) * sin(psi) * eta ^ 2 * q2 + 0.2e1 * cos(psi) * sin(psi) * q1 ^ 2 * q2 + cos(psi) ^ 2 * eta * q1 + cos(psi) * sin(psi) * eta * q2 - cos(psi) ^ 2 * q1 + 0.2e1 * cos(psi) * eta ^ 2 + 0.2e1 * cos(psi) * q1 ^ 2 + 0.2e1 * sin(psi) * q1 * q2 - q1 ^ 3 + 0.2e1 * eta * cos(psi) + eta * q1 + 0.2e1 * q1) / (0.1e1 + eta) / eta ^ 3;

% Angular Velocity in x direction with all the zonals and tesserals included
% sth = sin(f + g);
% cth = sin(f + g);
i = 2*atan(sqrt(p1^2 + p2^2));
% si = sin(i);
% ci = cos(i);
% b = a*sqrt(1 - e^2);
% mm = sqrt(mu/a^3);

w = AngVel([a,e,i,h,g,f]',theta,degree,order, mu, Re, Clm,Slm);

w__1 = w(1);

% RelPos with Psi matrix
RelPos = [-(q1 ^ 2 + q2 ^ 2 - 1) * sqrt(mu) / Vt * p ^ (-0.1e1 / 0.2e1) Vr * sqrt(mu) * sqrt(p) / Vt ^ 2 0 0 -a * (0.2e1 * sin(psi) * q1 * q2 + cos(psi) * (q1 ^ 2) - cos(psi) * (q2 ^ 2) + cos(psi) + (2 * q1)) / Vt ^ 2 / p * mu a * (sin(psi) * (q1 ^ 2) - sin(psi) * (q2 ^ 2) - 0.2e1 * cos(psi) * q1 * q2 - sin(psi) - (2 * q2)) / Vt ^ 2 / p * mu; 0 sqrt(p) / Vt * sqrt(mu) 0.2e1 * sqrt(p) * p2 / (p1 ^ 2 + p2 ^ 2 + 1) / Vt * sqrt(mu) -0.2e1 * sqrt(p) * p1 / (p1 ^ 2 + p2 ^ 2 + 1) / Vt * sqrt(mu) 0 0; 0 0 0.2e1 * sqrt(p) * sin(psi) / (p1 ^ 2 + p2 ^ 2 + 1) / Vt * sqrt(mu) -0.2e1 * sqrt(p) * cos(psi) / (p1 ^ 2 + p2 ^ 2 + 1) / Vt * sqrt(mu) 0 0;];

% RelVel with Psi Matrix
RelVel = [-Vr / a / 0.2e1 (2 * Vt) + sqrt(mu) * p ^ (-0.1e1 / 0.2e1) * (cos(psi) * p1 ^ 4 * q1 - 0.2e1 * cos(psi) * p1 ^ 2 * q1 * sigma1 ^ 2 - cos(psi) * p2 ^ 4 * q1 + sin(psi) * p1 ^ 4 * q2 - 0.2e1 * sin(psi) * p1 ^ 2 * q2 * sigma1 ^ 2 - sin(psi) * p2 ^ 4 * q2 - 0.2e1 * cos(psi) * q1 * sigma1 ^ 2 - 0.2e1 * sin(psi) * q2 * sigma1 ^ 2 + p1 ^ 4 - 0.2e1 * p1 ^ 2 * sigma1 ^ 2 - p2 ^ 4 - sigma1 ^ 4 - cos(psi) * q1 - sin(psi) * q2 - 0.4e1 * sigma1 ^ 2 - 0.2e1) / (sigma1 ^ 2 + 0.1e1) ^ 2 0 0 0.1e1 / p * a * q1 * Vr + sqrt(mu) * p ^ (-0.1e1 / 0.2e1) * sin(psi) 0.1e1 / p * a * q2 * Vr - sqrt(mu) * p ^ (-0.1e1 / 0.2e1) * cos(psi); -0.3e1 / 0.2e1 * Vt / a -Vr 0.2e1 * p2 / (sigma1 ^ 2 + 0.1e1) * Vr + 0.2e1 * w__1 * sqrt(p) * sin(psi) / (sigma1 ^ 2 + 0.1e1) / Vt * sqrt(mu) -0.2e1 * p1 / (sigma1 ^ 2 + 0.1e1) * Vr - 0.2e1 * w__1 * sqrt(p) * cos(psi) / (sigma1 ^ 2 + 0.1e1) / Vt * sqrt(mu) sqrt(mu) * p ^ (-0.3e1 / 0.2e1) * (cos(psi) * a * q1 ^ 2 + sin(psi) * a * q1 * q2 + cos(psi) * p + a * q1) + mu ^ (0.3e1 / 0.2e1) * a * (cos(psi) ^ 3 * q1 ^ 4 - 0.6e1 * cos(psi) ^ 3 * q1 ^ 2 * q2 ^ 2 + cos(psi) ^ 3 * q2 ^ 4 + 0.4e1 * cos(psi) ^ 2 * sin(psi) * q1 ^ 3 * q2 - 0.4e1 * cos(psi) ^ 2 * sin(psi) * q1 * q2 ^ 3 + cos(psi) ^ 3 * q1 ^ 2 - cos(psi) ^ 3 * q2 ^ 2 + 0.2e1 * cos(psi) ^ 2 * sin(psi) * q2 * q1 + 0.4e1 * cos(psi) ^ 2 * q1 ^ 3 - 0.8e1 * cos(psi) ^ 2 * q1 * q2 ^ 2 + 0.10e2 * cos(psi) * sin(psi) * q1 ^ 2 * q2 - 0.2e1 * cos(psi) * sin(psi) * q2 ^ 3 + 0.5e1 * q2 ^ 2 * cos(psi) * q1 ^ 2 - cos(psi) * q2 ^ 4 + 0.2e1 * sin(psi) * q1 * q2 ^ 3 + 0.2e1 * cos(psi) ^ 2 * q1 + 0.2e1 * cos(psi) * sin(psi) * q2 + 0.5e1 * cos(psi) * q1 ^ 2 + 0.6e1 * sin(psi) * q1 * q2 + 0.6e1 * q1 * q2 ^ 2 + cos(psi) + 0.2e1 * q1) * p ^ (-0.5e1 / 0.2e1) / (Vt ^ 2) sqrt(mu) * p ^ (-0.3e1 / 0.2e1) * (cos(psi) * a * q1 * q2 + sin(psi) * a * q2 ^ 2 + sin(psi) * p + a * q2) + mu ^ (0.3e1 / 0.2e1) * a * (-cos(psi) ^ 2 * sin(psi) * q2 ^ 4 - cos(psi) ^ 2 * sin(psi) * q2 ^ 2 + sin(psi) + 0.4e1 * q2 - 0.4e1 * q2 ^ 3 * cos(psi) ^ 2 - 0.2e1 * cos(psi) ^ 3 * q1 * q2 + 0.2e1 * sin(psi) * cos(psi) * q1 + 0.8e1 * cos(psi) * q1 * q2 - 0.2e1 * cos(psi) ^ 2 * q2 - sin(psi) * q1 ^ 2 + 0.6e1 * sin(psi) * q2 ^ 2 + sin(psi) * q2 ^ 4 + 0.4e1 * q2 ^ 3 + sin(psi) * cos(psi) ^ 2 * q1 ^ 2 - sin(psi) * q1 ^ 2 * q2 ^ 2 + 0.4e1 * cos(psi) * q1 * q2 ^ 3 + 0.6e1 * cos(psi) ^ 2 * sin(psi) * q1 ^ 2 * q2 ^ 2 + 0.10e2 * cos(psi) * sin(psi) * q1 * q2 ^ 2 + 0.4e1 * cos(psi) ^ 3 * q1 ^ 3 * q2 - 0.4e1 * cos(psi) ^ 3 * q1 * q2 ^ 3 - cos(psi) ^ 2 * sin(psi) * q1 ^ 4 + 0.8e1 * cos(psi) ^ 2 * q1 ^ 2 * q2 - 0.2e1 * cos(psi) * sin(psi) * q1 ^ 3 - 0.2e1 * q2 * cos(psi) * q1 ^ 3 - 0.2e1 * q1 ^ 2 * q2) * p ^ (-0.5e1 / 0.2e1) / (Vt ^ 2); 0 -w__1 * sqrt(p) / Vt * sqrt(mu) 0.2e1 * sin(psi) / (sigma1 ^ 2 + 0.1e1) * Vr + 0.2e1 / (sigma1 ^ 2 + 0.1e1) * cos(psi) * Vt - 0.2e1 * w__1 * sqrt(p) * p2 / (sigma1 ^ 2 + 0.1e1) / Vt * sqrt(mu) -0.2e1 / (sigma1 ^ 2 + 0.1e1) * cos(psi) * Vr - 0.2e1 * sqrt(mu) * p ^ (-0.1e1 / 0.2e1) * (cos(psi) ^ 2 * q2 - sin(psi) * cos(psi) * q1 - sin(psi) - q2) / (sigma1 ^ 2 + 0.1e1) + 0.2e1 * w__1 * sqrt(p) * p1 / (sigma1 ^ 2 + 0.1e1) / Vt * sqrt(mu) 0 0;];

if DiffTA == true

% Delta RelPos Matrix for Lambda
DelRelPosm = [0 dTLML * Vr * sqrt(mu) * sqrt(p) / Vt ^ 2 - Vr * sqrt(mu) * sqrt(p) / Vt ^ 2 0 0 Vr * sqrt(mu) * sqrt(p) / Vt ^ 2 * dTLq1 Vr * sqrt(mu) * sqrt(p) / Vt ^ 2 * dTLq2; 0 dTLML * sqrt(p) / Vt * sqrt(mu) - sqrt(p) / Vt * sqrt(mu) 0 0 sqrt(p) / Vt * sqrt(mu) * dTLq1 sqrt(p) / Vt * sqrt(mu) * dTLq2; 0 0 0 0 0 0;];

% Delta RelVel Matrix for Lambda
DelRelVelm = [0 dTLML * ((2 * Vt) + sqrt(mu) * p ^ (-0.1e1 / 0.2e1) * (cos(psi) * p1 ^ 4 * q1 - 0.2e1 * cos(psi) * p1 ^ 2 * q1 * sigma1 ^ 2 - cos(psi) * p2 ^ 4 * q1 + sin(psi) * p1 ^ 4 * q2 - 0.2e1 * sin(psi) * p1 ^ 2 * q2 * sigma1 ^ 2 - sin(psi) * p2 ^ 4 * q2 - 0.2e1 * cos(psi) * q1 * sigma1 ^ 2 - 0.2e1 * sin(psi) * q2 * sigma1 ^ 2 + p1 ^ 4 - 0.2e1 * p1 ^ 2 * sigma1 ^ 2 - p2 ^ 4 - sigma1 ^ 4 - cos(psi) * q1 - sin(psi) * q2 - 0.4e1 * sigma1 ^ 2 - 0.2e1) / (sigma1 ^ 2 + 0.1e1) ^ 2) - (2 * Vt) - sqrt(mu) * p ^ (-0.1e1 / 0.2e1) * (cos(psi) * p1 ^ 4 * q1 - 0.2e1 * cos(psi) * p1 ^ 2 * q1 * sigma1 ^ 2 - cos(psi) * p2 ^ 4 * q1 + sin(psi) * p1 ^ 4 * q2 - 0.2e1 * sin(psi) * p1 ^ 2 * q2 * sigma1 ^ 2 - sin(psi) * p2 ^ 4 * q2 - 0.2e1 * cos(psi) * q1 * sigma1 ^ 2 - 0.2e1 * sin(psi) * q2 * sigma1 ^ 2 + p1 ^ 4 - 0.2e1 * p1 ^ 2 * sigma1 ^ 2 - p2 ^ 4 - sigma1 ^ 4 - cos(psi) * q1 - sin(psi) * q2 - 0.4e1 * sigma1 ^ 2 - 0.2e1) / (sigma1 ^ 2 + 0.1e1) ^ 2 0 0 ((2 * Vt) + sqrt(mu) * p ^ (-0.1e1 / 0.2e1) * (cos(psi) * p1 ^ 4 * q1 - 0.2e1 * cos(psi) * p1 ^ 2 * q1 * sigma1 ^ 2 - cos(psi) * p2 ^ 4 * q1 + sin(psi) * p1 ^ 4 * q2 - 0.2e1 * sin(psi) * p1 ^ 2 * q2 * sigma1 ^ 2 - sin(psi) * p2 ^ 4 * q2 - 0.2e1 * cos(psi) * q1 * sigma1 ^ 2 - 0.2e1 * sin(psi) * q2 * sigma1 ^ 2 + p1 ^ 4 - 0.2e1 * p1 ^ 2 * sigma1 ^ 2 - p2 ^ 4 - sigma1 ^ 4 - cos(psi) * q1 - sin(psi) * q2 - 0.4e1 * sigma1 ^ 2 - 0.2e1) / (sigma1 ^ 2 + 0.1e1) ^ 2) * dTLq1 ((2 * Vt) + sqrt(mu) * p ^ (-0.1e1 / 0.2e1) * (cos(psi) * p1 ^ 4 * q1 - 0.2e1 * cos(psi) * p1 ^ 2 * q1 * sigma1 ^ 2 - cos(psi) * p2 ^ 4 * q1 + sin(psi) * p1 ^ 4 * q2 - 0.2e1 * sin(psi) * p1 ^ 2 * q2 * sigma1 ^ 2 - sin(psi) * p2 ^ 4 * q2 - 0.2e1 * cos(psi) * q1 * sigma1 ^ 2 - 0.2e1 * sin(psi) * q2 * sigma1 ^ 2 + p1 ^ 4 - 0.2e1 * p1 ^ 2 * sigma1 ^ 2 - p2 ^ 4 - sigma1 ^ 4 - cos(psi) * q1 - sin(psi) * q2 - 0.4e1 * sigma1 ^ 2 - 0.2e1) / (sigma1 ^ 2 + 0.1e1) ^ 2) * dTLq2; 0 -dTLML * Vr + Vr 0 0 -Vr * dTLq1 -Vr * dTLq2; 0 -dTLML * w__1 * sqrt(p) / Vt * sqrt(mu) + w__1 * sqrt(p) / Vt * sqrt(mu) 0 0 -w__1 * sqrt(p) / Vt * sqrt(mu) * dTLq1 -w__1 * sqrt(p) / Vt * sqrt(mu) * dTLq2;];

else
    DelRelPosm = 0;
    DelRelVelm = 0;
end

Sigma = [RelPos + DelRelPosm;
         RelVel + DelRelVelm];

end

