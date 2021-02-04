% Mean2TrueLong: Mean to True Longitude Solver using Newton-Raphson method
%   Note: MA in radian, MAX Iterations: 1000
% Author: Bharat Mahajan (https://github.com/princemahajan)

function [TL, EL] = Mean2TrueLong(ML, q1,q2, tol)

% Yes 1000 not 999!
MAX_ITERATIONS = 1000;

% initial guess
EL = ML;

ctr = 0;

f = 9999;
dEL = 9999;

% Newton_Raphson
while abs(dEL) > tol && abs(f) > tol && (ctr < MAX_ITERATIONS)
    
    [f, dfdx] = ModKeplerEqn(EL,ML,q1,q2);
    
    dEL = -f / dfdx;
    
    EL = EL + dEL;
    ctr = ctr + 1;
end

if ctr == MAX_ITERATIONS
    disp('Mean to True Longitude Solver failed...');
    EL = NaN;
    return;
end

% Make sure EL is between 0 and 2*pi

EL = mod(mod(EL,2*pi),2*pi);
  
% calculate true longitude
abyab = 1/(1 + sqrt(1 - q1^2 - q2^2));

abyr = 1/(1 - q1*cos(EL) - q2*sin(EL));

sP = abyr*((1 - abyab*q1^2)*sin(EL) + abyab*q1*q2*cos(EL) - q2);
cP = abyr*((1 - abyab*q2^2)*cos(EL) + abyab*q1*q2*sin(EL) - q1);

TL = atan2(sP, cP);

TL = mod(mod(TL,2*pi),2*pi);

end

function [f, dfdEL] = ModKeplerEqn(EL,ML,q1,q2)

f = ML - (EL + q2*cos(EL) - q1*sin(EL));

dfdEL = - (1 - q1*cos(EL) - q2*sin(EL));

end
