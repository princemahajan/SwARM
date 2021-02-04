%% PCO Formation Differential Initial Elements
%   Desc: Computes equinoctial initial conditions for Deputy. 
% Author: Bharat Mahajan (https://github.com/princemahajan)

function DelEq = GCODiffEquinoctial(rhox,alphax1,rhoy,rhoz,alphaz1,ChiefEq)

a = ChiefEq(1);
ML = ChiefEq(2);
p1 = ChiefEq(3);
p2 = ChiefEq(4);
q1 = ChiefEq(5);
q2 = ChiefEq(6);

RAAN = atan2(p2,p1);

% adjust phase angles
alphax = alphax1 - RAAN;
alphaz = alphaz1 - RAAN;

% p1,p2

dp1 = (1 + p1^2 + p2^2)/2/a*rhoz*cos(alphaz);
dp2 = -(1 + p1^2 + p2^2)/2/a*rhoz*sin(alphaz);

% q1,q2
dq1 = -rhox/a*sin(alphax);
dq2 = -rhox/a*cos(alphax);

% ML
dML = rhoy/a - 2*(p2*dp1 - p1*dp2)/(1 + p1^2 + p2^2);

dML1 = -rhoy/(a*(1+2*(sin(ML)*q2 + cos(ML)*q1)));


DelEq = [0,dML,dp1,dp2,dq1,dq2]';

end
