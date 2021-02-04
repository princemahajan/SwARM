% KeplerEqSolver: Kepler Equation Solver using Newton-Raphson method
%   Note: MA in radian, MAX Iterations: 999
% Author: Bharat Mahajan (https://github.com/princemahajan)

function [EA, TA] = KeplerEqSolver(MA, ecc, tol)

% Yes 1000 not 999!
MAX_ITERATIONS = 1000;

% compute the correct quadrant of MA
FullRevs = fix(MA/2/pi);

% primary MA value
MA1 = mod(mod(MA,2*pi),2*pi);
if abs(MA1) < tol
    MA1 = 0;
end

% initial guess
EA1 = MA1 + (ecc/2); % see Prussing and Conway book

ctr = 0;

deltaM = 9999;

% Newton_Raphson
while (abs(deltaM) > tol) && (ctr <= MAX_ITERATIONS)
    dMdE = 1 - ecc * cos(EA1);
    deltaM = (EA1 - ecc * sin(EA1)) - MA1;
    deltaE = -deltaM / dMdE;
    EA1 = EA1 + deltaE;
    ctr = ctr + 1;
end

if ctr == MAX_ITERATIONS
    disp('KeplerEqSolver failed...');
    EA = NaN; TA = NaN;
    return;
end

% calculate true anomaly using half angle formula
eterm = sqrt((1+ecc)/(1-ecc));
TA1 = 2*atan(eterm * tan(EA1/2));

% make sure TA,EA is in [0,2*pi), yes 2 times mod!
EA1 = mod(mod(EA1,2*pi),2*pi);
TA1 = mod(mod(TA1,2*pi),2*pi);

% Use exact zero for below tolerance values
if abs(TA1) < tol
    TA1 = 0;
end
if abs(EA1) < tol
    EA1 = 0;
end

% adjust for negative anomalies
if TA1 ~= 0 && MA < 0
    TA1 = -2*pi + TA1;
end
if EA1 ~= 0 && MA < 0
    EA1 = -2*pi + EA1;
end

% anomalies in the correct quadrant
TA = TA1 + FullRevs*2*pi;
EA = EA1 + FullRevs*2*pi;

end
