%% EqLPZonal
%     Desc: Computes long-period effects up to first order for Jn (n > 2) for Equinoctial elements.
%     Equinoctial elements: [a l+g+h tan(i/2)*cos(h) tan(i/2)*sin(h) e*cos(g+h) e*sin(g+h)];
% Reference: 
% B. Mahajan, S. R. Vadali and K. T. Alfriend, �Analytic Solution For Satellite 
% Relative Motion: The Complete Zonal Gravitational Problem,� 26th AAS/AIAA 
% Space Flight Mechanics Meeting, Napa, CA, February 2016.
% Author: Bharat Mahajan (https://github.com/princemahajan)

function [Eqlp, Dlp] = EqLPZonal(Xm,degree,mu,Re,Jcoeff,InverseOn, JacobianOn,tol) 

O2Msign = 1;
if InverseOn == true
	O2Msign = -1;
end

% for each zonal
narr = 3:1:degree;
J2 = Jcoeff(2);

a = Xm(1); Lambda = Xm(2); p1 = Xm(3); p2 = Xm(4); q1 = Xm(5); q2 = Xm(6);

DelXm = Del2Eqn(Xm,mu,true);
 
l = DelXm(1); g = DelXm(2); h = DelXm(3);L = DelXm(4);G =DelXm(5); H= DelXm(6);

e = sqrt(q1^2 + q2^2);
eta = sqrt(1 - e^2);
i = 2*atan(sqrt(p1^2 + p2^2));
si = sin(i); 
ci = cos(i);
ch = cos(h); 
sh = sin(h);
si2 = sin(i/2);
ci2 = cos(i/2);

sgh = sin(g+h); cgh = cos(g+h);
sg = sin(g); cg = cos(g);

% J2 Equinoctial First Order Long-Period Effects (Using Maple)
J = J2;
R__e = Re;

% First order Long-Period effects due to higher degree zonals
Jn_LPX1 = zeros(6,1);
Dlpn = zeros(6);

for ctr = 1:length(narr)
    
    n = narr(ctr);
    Jn = Jcoeff(n);
    
    % SMA
    anlp1 = 0;
    
    % Mean Longitude
    [lW1W2] = getW1W2(@MlongTerm,n,Jn,J2,a,e,i,eta,g,mu,Re);
    Lambdanlp1 = sin(i)/sqrt(mu*a)*lW1W2;
    
    % p1,p2
    [pW1W2] = getW1W2(@p1p2Term,n,Jn,J2,a,e,i,eta,g,mu,Re);
    [pW1W4] = getW1W4(@dummy,n,Jn,J2,a,e,i,eta,g,mu,Re);
    p1nlp1 = cos(i)/(2*G*(cos(i/2))^2)*(sin(h)*pW1W2 - cos(h)*pW1W4);
    p2nlp1 = cos(i)/(2*G*(cos(i/2))^2)*(-cos(h)*pW1W2 - sin(h)*pW1W4);
    
    % q1,q2
    [qW1W2] = getW1W2(@q1q2Term1,n,Jn,J2,a,e,i,eta,g,mu,Re);
    [qW3W2] = getW3W2(@q1q2Term2,n,Jn,J2,a,e,i,eta,g,mu,Re);
    [qW3W4] = getW3W4(@dummy,n,Jn,J2,a,e,i,eta,g,mu,Re);
    q1nlp1 = sin(i)*sin(g+h)/G*(qW1W2 + eta^2*qW3W2) + sin(i)*eta/L*cos(g+h)*qW3W4;
    q2nlp1 = -sin(i)*cos(g+h)/G*(qW1W2 + eta^2*qW3W2) + sin(i)*eta/L*sin(g+h)*qW3W4;
    
    Jn_LPX1 = Jn_LPX1 + J2*[anlp1; Lambdanlp1; p1nlp1; p2nlp1; q1nlp1; q2nlp1];
    
    
    if JacobianOn ~= true
        continue;
    end
    
    % SMA LP partials
    dadx = [0,0,0,0,0,0];
    
    % Mean Longitude partials
    dMLda = (-1/2/a)*Lambdanlp1 -(2*n-5)/2/a*Lambdanlp1;
    
    dMLdML = 0;
    
    W1W2t23 = getW1W2(@MLp12termW12,n,Jn,J2,a,e,i,eta,g,mu,Re);
    W1W4t4 = getW1W4(@MLp12termW14,n,Jn,J2,a,e,i,eta,g,mu,Re);

    dMLdp1 = 1/sqrt(mu*a)*(cos(i)*lW1W2*(1+ci)*ch + W1W2t23*(1+ci)*ch + W1W4t4*(1+ci)*sh); 
    dMLdp2 = 1/sqrt(mu*a)*(cos(i)*lW1W2*(1+ci)*sh + W1W2t23*(1+ci)*sh - W1W4t4*(1+ci)*ch); 
    
    mldq12t1 = getW3W2(@MLq12termW321,n,Jn,J2,a,e,i,eta,g,mu,Re);
    mldq12t2 = getW3W2(@MLq12termW322,n,Jn,J2,a,e,i,eta,g,mu,Re);
    mldq12t3 = getW3W4(@MLq12termW343,n,Jn,J2,a,e,i,eta,g,mu,Re);
    
    dMLdq1 = si/sqrt(mu*a)*(mldq12t1*q1 + mldq12t2*cgh + mldq12t3*(-sgh));
    dMLdq2 = si/sqrt(mu*a)*(mldq12t1*q2 + mldq12t2*sgh + mldq12t3*(cgh));
    
    % p1
    dp1da = -(2*n - 4)/2/a*p1nlp1;
    
    dp1dML = 0;
    
    t32 = e*getW3W2(@p12p12term32,n,Jn,J2,a,e,i,eta,g,mu,Re);
    t34 = e*getW3W4(@p12p12term34,n,Jn,J2,a,e,i,eta,g,mu,Re);
    t3b4 = e*getW3bW4(@p12p12term3b4,n,Jn,J2,a,e,i,eta,g,mu,Re);
    t3b2 = e*getW3bW2(@p12p12term3b2,n,Jn,J2,a,e,i,eta,g,mu,Re);
    t3bnjp = e*getW3bW2njp(@p12p12term3b2njp,n,Jn,J2,a,e,i,eta,g,mu,Re);
    t3bjstar = e*getW3jstar(@p12p12termjstar,(n-1)/2,n,Jn,J2,a,e,i,eta,g,mu,Re);
    
    dp1dp1 = 1/eta/sqrt(mu*a)*(sh*t32*ch + ch*t34*ch + sh*t3b4*sh + ch*t3b2*sh + ch*t3bnjp*sh + (sh*sg+ch*cg)*t3bjstar*sh);
    
    dp1dp2 = 1/eta/sqrt(mu*a)*(sh*t32*sh + ch*t34*sh - sh*t3b4*ch - ch*t3b2*ch - ch*t3bnjp*ch - (sh*sg+ch*cg)*t3bjstar*ch);

    t32sc = getW3W2(@p1q1Term32shcgh,n,Jn,J2,a,e,i,eta,g,mu,Re);
    t34cc = getW3W4(@p1q1Term34chcgh,n,Jn,J2,a,e,i,eta,g,mu,Re);
    t34ss = getW3W4(@p1q1Term34shsgh,n,Jn,J2,a,e,i,eta,g,mu,Re);
    t32njpcs = getW3W2njp(@p1q1Term32njpchsgh,n,Jn,J2,a,e,i,eta,g,mu,Re);
    
    dp1dq1 = 1/eta^2*p1nlp1*q1 + ci/eta/sqrt(mu*a)/(1+ci)*(sh*t32sc*cgh + ch*t34cc*cgh + sh*t34ss*sgh + ch*t32njpcs*sgh);
    
    dp1dq2 = 1/eta^2*p1nlp1*q2 + ci/eta/sqrt(mu*a)/(1+ci)*(sh*t32sc*sgh + ch*t34cc*sgh - sh*t34ss*cgh - ch*t32njpcs*cgh);
    
    % p2
    dp2da = -(2*n - 4)/2/a*p2nlp1;
    
    dp2dML = 0;
    
    dp2dp1 = 1/eta/sqrt(mu*a)*(-ch*t32*ch + sh*t34*ch - ch*t3b4*sh + sh*t3b2*sh + sh*t3bnjp*sh - (ch*sg-sh*cg)*t3bjstar*sh);
    
    dp2dp2 = 1/eta/sqrt(mu*a)*(-ch*t32*sh + sh*t34*sh + ch*t3b4*ch - sh*t3b2*ch - sh*t3bnjp*ch + (ch*sg-sh*cg)*t3bjstar*ch);
    
    dp2dq1 = 1/eta^2*p2nlp1*q1 + ci/eta/sqrt(mu*a)/(1+ci)*(-ch*t32sc*cgh + sh*t34cc*cgh - ch*t34ss*sgh + sh*t32njpcs*sgh);
    
    dp2dq2 = 1/eta^2*p2nlp1*q2 + ci/eta/sqrt(mu*a)/(1+ci)*(-ch*t32sc*sgh + sh*t34cc*sgh + ch*t34ss*cgh - sh*t32njpcs*cgh);
    
    % q1
    dq1da = -(2*n-4)/2/a*q1nlp1;
    
    dq1dML = 0;
    
    t12sghch = getW3W2(@q1p1Term32sghch,n,Jn,J2,a,e,i,eta,g,mu,Re);
    t14cghch = getW3W4(@q1p1Term34cghch,n,Jn,J2,a,e,i,eta,g,mu,Re);
    t14sghsh = getW3W4(@q1p1Term34sghsh,n,Jn,J2,a,e,i,eta,g,mu,Re);
    t32njpcghsh = getW3W2(@q1p1Term32njpcghsh,n,Jn,J2,a,e,i,eta,g,mu,Re);
    
    dq1dp1 = 1/eta/sqrt(mu*a)*(t12sghch*sgh*ch + eta^2*t14cghch*cgh*ch + t14sghsh*sgh*sh  ...
                    + eta^2*t32njpcghsh*cgh*sh); 
    
    dq1dp2 = 1/eta/sqrt(mu*a)*(t12sghch*sgh*sh + eta^2*t14cghch*cgh*sh - t14sghsh*sgh*ch  ...
                    - eta^2*t32njpcghsh*cgh*ch); 
                
    t32sc = getW3W2(@q1q1Term32sghcgh,n,Jn,J2,a,e,i,eta,g,mu,Re);
    t322sc = getW3bW2(@q1q1Term322sghcgh,n,Jn,J2,a,e,i,eta,g,mu,Re);
    t34cc = getW3W4(@q1q1Term34cghcgh,n,Jn,J2,a,e,i,eta,g,mu,Re);
    t324cc = getW3bW4(@q1q1Term324cghcgh,n,Jn,J2,a,e,i,eta,g,mu,Re);
    t32cs = getW3W2(@q1q1Term32cghsgh,n,Jn,J2,a,e,i,eta,g,mu,Re);
    t34ss = getW3W4(@q1q1Term34sghsgh,n,Jn,J2,a,e,i,eta,g,mu,Re);
    t3b4ss = getW3bW4(@q1q1Term3b4sghsgh,n,Jn,J2,a,e,i,eta,g,mu,Re);
    t3b2cs = getW3bW2(@q1q1Term3b2cghsgh,n,Jn,J2,a,e,i,eta,g,mu,Re);
    t3b2njpcs = getW3bW2njp(@q1q1Term3b2njpcghsgh,n,Jn,J2,a,e,i,eta,g,mu,Re);
           
    dq1dq1 = 1/eta/sqrt(mu*a)*(sgh*t32sc*cgh + sgh*t322sc*cgh + cgh*t34cc*cgh + cgh*t324cc*cgh ...
                + cgh*t32cs*sgh + sgh*t34ss*sgh + sgh*t3b4ss*sgh + cgh*t3b2cs*sgh + cgh*t3b2njpcs*sgh);
    
    dq1dq2 = 1/eta/sqrt(mu*a)*(sgh*t32sc*sgh + sgh*t322sc*sgh + cgh*t34cc*sgh + cgh*t324cc*sgh ...
                - cgh*t32cs*cgh - sgh*t34ss*cgh - sgh*t3b4ss*cgh - cgh*t3b2cs*cgh - cgh*t3b2njpcs*cgh);
    
    % q2
    
    dq2da = -(2*n-4)/2/a*q2nlp1;
    
    dq2dML = 0;
    
    dq2dp1 = 1/eta/sqrt(mu*a)*(-t12sghch*cgh*ch + eta^2*t14cghch*sgh*ch - t14sghsh*cgh*sh  ...
                    + eta^2*t32njpcghsh*sgh*sh);
                
    dq2dp2 = 1/eta/sqrt(mu*a)*(-t12sghch*cgh*sh + eta^2*t14cghch*sgh*sh + t14sghsh*cgh*ch  ...
                    - eta^2*t32njpcghsh*sgh*ch);

    dq2dq1 = 1/eta/sqrt(mu*a)*(-cgh*t32sc*cgh -cgh*t322sc*cgh + sgh*t34cc*cgh + sgh*t324cc*cgh ...
                + sgh*t32cs*sgh - cgh*t34ss*sgh - cgh*t3b4ss*sgh + sgh*t3b2cs*sgh + sgh*t3b2njpcs*sgh);

    dq2dq2 = 1/eta/sqrt(mu*a)*(-cgh*t32sc*sgh -cgh*t322sc*sgh + sgh*t34cc*sgh + sgh*t324cc*sgh ...
                - sgh*t32cs*cgh + cgh*t34ss*cgh + cgh*t3b4ss*cgh - sgh*t3b2cs*cgh - sgh*t3b2njpcs*cgh);


    DLPn = [dadx;
            dMLda, dMLdML, dMLdp1, dMLdp2, dMLdq1, dMLdq2;
            dp1da, dp1dML, dp1dp1, dp1dp2, dp1dq1, dp1dq2;
            dp2da, dp2dML, dp2dp1, dp2dp2, dp2dq1, dp2dq2;
            dq1da, dq1dML, dq1dp1, dq1dp2, dq1dq1, dq1dq2;
            dq2da, dq2dML, dq2dp1, dq2dp2, dq2dq1, dq2dq2];            
            
    Dlpn = Dlpn + DLPn;
        
end

% Long-period effects
EqLPDelX1 = O2Msign*Jn_LPX1;
EqLPDelX2 = O2Msign*0;
Eqlp = EqLPDelX1 + EqLPDelX2;

% Dlp matrix
Dlp = O2Msign*J2*Dlpn;

end



% Intermediate term for Mean Longitude
function func = MlongTerm(k,j,n,i,e,eta)

func = -(2*n-5)/eta -eta*k/(1+eta) + 10*cos(i)*(1 - cos(i))/(eta*(1 - 5*(cos(i))^2)) ...
            + (n-2*j)/eta*(-cos(i)/2/(cos(i/2))^2);

end

% Intermediate term for Mean Longitude partials wrt p1p2 term2-term3
function func = MLp12termW12(k,j,n,i,e,eta)

ci = cos(i); si=sin(i);

func = (-(2*n-5)/eta -eta*k/(1+eta) + 10*ci*(1-ci)/(eta*(1-5*ci^2)) ...
            + (n-2*j)/eta*(-ci/(1+ci)))*(-10*si^2*ci/(1-5*ci^2) + (n-2*j-1)*ci) ...
            + si^2/eta/(1-5*ci^2)*((-10+20*ci) - 100*ci^2*(1-ci)/(1-5*ci^2)) ...
            + si/eta/(1+ci)*(si - ci*tan(i/2));

end

function func = MLp12termW14(k,j,n,i,e,eta)

ci = cos(i); 

func = (-(2*n-5)/eta -eta*k/(1+eta) + 10*ci*(1-ci)/(eta*(1-5*ci^2)) ...
            + (n-2*j)/eta*(-ci/(1+ci)));

end

% Intermediate term for Mean Longitude partials wrt q1q2 term2-term3
function func = MLq12termW321(k,j,n,i,e,eta)

ci = cos(i); si=sin(i);

func = (-(2*n-5)/eta -eta*k/(1+eta) + 10*ci*(1-ci)/(eta*(1-5*ci^2)) ...
            + (n-2*j)/eta*(-ci/(1+ci)))*(2*n-5)*e/eta^2 + e*(-(2*n-5)/eta^3 +k/eta/(1+eta)^2 ...
            + 10*ci*(1-ci)/(eta^3*(1-5*ci^2) + (n-2*j)/eta^3*(-ci/(1+ci))));

end

function func = MLq12termW322(k,j,n,i,e,eta)

ci = cos(i); si=sin(i);

func = (-(2*n-5)/eta -eta*k/(1+eta) + 10*ci*(1-ci)/(eta*(1-5*ci^2)) ...
            + (n-2*j)/eta*(-ci/(1+ci)))*k;

end

function func = MLq12termW343(k,j,n,i,e,eta)

ci = cos(i); si=sin(i);

func = (-(2*n-5)/eta -eta*k/(1+eta) + 10*ci*(1-ci)/(eta*(1-5*ci^2)) ...
            + (n-2*j)/eta*(-ci/(1+ci)));

end

% Intermediate terms for p12 partials wrt p12



function term = p12p12term32(k,j,n,i,e,eta)

si = sin(i);ci = cos(i);

term = tan(i/2)*((n-2*j - 10*(si)^2/(1 - 5*ci^2))*(-1 -10*ci^2/(1-5*ci^2)*(1+ci)) ...
        - (n-2*j-1)*10*ci^2/(1-5*ci^2)*(1+ci) + 80*ci^2/(1-5*ci^2)^2*(1+ci))...
            + (n-2*j)*(n-2*j-1)*ci^2;

end

function term = p12p12term34(k,j,n,i,e,eta)

si = sin(i);ci = cos(i);

term = tan(i/2)*(1 + 10*ci^2/(1-5*ci^2)*(1+ci)) - (n-2*j-1)*ci^2;

end

function term = p12p12term3b4(k,j,n,i,e,eta)

si = sin(i);ci = cos(i);

term = (-1 + (n-2*j - 10*(si)^2/(1 - 5*ci^2)))*ci;

end

function term = p12p12term3b2(k,j,n,i,e,eta)

si = sin(i);ci = cos(i);

term = (-(n-2*j - 10*(si)^2/(1 - 5*ci^2)))*ci;

end

function term = p12p12term3b2njp(k,j,n,i,e,eta)

si = sin(i);ci = cos(i);

term = 1*ci;

end

function term = p12p12termjstar(k,j,n,i,e,eta)

if mod(n,2) == 0 
    term = 0;
else

si = sin(i);ci = cos(i);

term = ci*10*si/(1-5*ci^2)*getGamma(0,j,n)*nchoosek(k, (k-1)/2);
end

end

function term = p1q1Term32shcgh(k,j,n,i,e,eta)

si = sin(i);ci = cos(i);

term = (n-2*j - 10*(si)^2/(1 - 5*ci^2))*((2*n-5)*e^2/eta^2 + k);

end

function term = p1q1Term34chcgh(k,j,n,i,e,eta)

si = sin(i);ci = cos(i);

term = (-(2*n-5)*e^2/eta^2 - k);

end

function term = p1q1Term34shsgh(k,j,n,i,e,eta)

si = sin(i);ci = cos(i);

term = -(n-2*j - 10*(si)^2/(1 - 5*ci^2));

end

function term = p1q1Term32njpchsgh(k,j,n,i,e,eta)

si = sin(i);ci = cos(i);

term = -1;

end

% intermediate terms for q1 q2 partials

function term = q1p1Term32sghch(k,j,n,i,e,eta)

si = sin(i); ci = cos(i);

termi = (e*(2*n-5) - 5*e*sin(2*i)*tan(i/2)/(1-5*ci^2) + e*ci*(n-2*j)/(1+ci));

term = termi*(ci*(1+ci)*e + e*(-10*si^2*ci/(1 - 5*ci^2) + (n-2*j-1)*ci)*(1+ci)) ...
        + e*(10*si^2*e*(5*ci^2 -2*ci +1)/(1-5*ci^2)^2 + e*si*ci*(n-2*j)*tan(i/2)/(1+ci))*(1+ci) ...
        + eta^2*ci*k*(1+ci) + (-10*si^2*ci/(1 - 5*ci^2) + (n-2*j-1)*ci)*k*(1+ci);

end

function term = q1p1Term34cghch(k,j,n,i,e,eta)

si = sin(i); ci = cos(i);

term = ci*(1+ci) + (1+ci)*(-10*si^2*ci/(1 - 5*ci^2) + (n-2*j-1)*ci);

end

function term = q1p1Term34sghsh(k,j,n,i,e,eta)

si = sin(i); ci = cos(i);

termi = (e*(2*n-5) - 5*e*sin(2*i)*tan(i/2)/(1-5*ci^2) + e*ci*(n-2*j)/(1+ci));

term = e*termi*(1+ci) + eta^2*k*(1+ci);

end

function term = q1p1Term32njpcghsh(k,j,n,i,e,eta)

si = sin(i); ci = cos(i);

term = -(1+ci);

end

function term = q1q1Term32sghcgh(k,j,n,i,e,eta)

si = sin(i); ci = cos(i);

t1 = ((2*n-5) - 5*sin(2*i)*tan(i/2)/(1-5*ci^2) + ci*(n-2*j)/(1+ci));

term = t1*(si*e^2/eta^2*e + (2*n-5)*e^2/eta^2*si + e*si + e*si) + si*e/eta^2*k -2*e*k*si + (2*n-5)*e*si*k;

end

function term = q1q1Term322sghcgh(k,j,n,i,e,eta)

si = sin(i); ci = cos(i);

t1 = ((2*n-5) - 5*sin(2*i)*tan(i/2)/(1-5*ci^2) + ci*(n-2*j)/(1+ci));

term = (t1*(e*si) + (eta^2*k*si))*(k-1);

end


function term = q1q1Term34cghcgh(k,j,n,i,e,eta)

si = sin(i); ci = cos(i);

term = si*(-e) + si*(2*n-5)*e;

end

function term = q1q1Term324cghcgh(k,j,n,i,e,eta)

si = sin(i); ci = cos(i);

term = si*eta^2*(k-1);

end

function term = q1q1Term32cghsgh(k,j,n,i,e,eta)

si = sin(i); ci = cos(i);

t1 = e*((2*n-5) - 5*sin(2*i)*tan(i/2)/(1-5*ci^2) + ci*(n-2*j)/(1+ci));

term = -t1*si;

end
function term = q1q1Term34sghsgh(k,j,n,i,e,eta)

si = sin(i); ci = cos(i);
t1 = e*((2*n-5) - 5*sin(2*i)*tan(i/2)/(1-5*ci^2) + ci*(n-2*j)/(1+ci));
term = -si*e*t1;

end
function term = q1q1Term3b4sghsgh(k,j,n,i,e,eta)

si = sin(i); ci = cos(i);

term = -si*(eta^2*k - eta^2);

end

function term = q1q1Term3b2cghsgh(k,j,n,i,e,eta)

si = sin(i); ci = cos(i);

term = -si*(eta^2*k);

end

function term = q1q1Term3b2njpcghsgh(k,j,n,i,e,eta)

si = sin(i); ci = cos(i);

term = -si*(-eta^2);

end

% Intermediate term for p1,p2
function func = p1p2Term(k,j,n,i,e,eta)

func = n-2*j - 10*(sin(i))^2/(1 - 5*(cos(i))^2);

end

% Intermediate term-1 for q1q2
function func = q1q2Term1(k,j,n,i,e,eta)

func = e*(2*n-5) - 5*e*sin(2*i)*tan(i/2)/(1-5*(cos(i))^2) + e*cos(i)*(n-2*j)/2/((cos(i/2))^2);

end

% Intermediate term-2 for q1q2
function func = q1q2Term2(k,j,n,i,e,eta)

func = k;

end

% dummy
function func = dummy(k,j,n,i,e,eta)

func = 1;

end


function [W1FW2] = getW1W2(func,n,Jn,J2,a,e,i,eta,g,mu,Re)

deln = deltalp(n,Jn,J2,a,eta,mu,Re);

sumj1 = 0;

% outer j loop
for j = 0:floor(n/2)
   
    blp = betalp(j,n,i);

    sumk1 = 0;
    

    % innner k loop
    for k1 = 0:(floor(n/2)-1)
        
        if mod(n,2) == 0
            % even harmonic
            k = 2*k1;
        else
            % odd harmonic
            k = 2*k1 + 1;
        end
        
        alp = getAlpha(k,n,e);
        
        % compute W2, W4
        [W2] = getW2(n,j,k,g);
        
        sumk1 = sumk1 + alp/2^k*func(k,j,n,i,e,eta)*W2;
        
    end
    
    if isinf(blp) == 0
        sumj1 = sumj1 + blp*sumk1;
    elseif sumk1 ~= 0
        disp('Inclination singularity');
    end
        
end

W1FW2 = deln*sumj1;

end

function [W1FW2] = getW1W2jk(func,jf,k0,n,Jn,J2,a,e,i,eta,g,mu,Re)

deln = deltalp(n,Jn,J2,a,eta,mu,Re);

sumj1 = 0;

% outer j loop
for j = 0:floor((n-jf)/2)
   
    blp = betalpjf(j,n,i,jf);

    sumk1 = 0;
    

    % innner k loop
    for k1 = (k0/2):(floor(n/2)-1)
        
        if mod(n,2) == 0
            % even harmonic
            k = 2*k1;
        else
            % odd harmonic
            k = 2*k1 + 1;
        end
        
        alp = getAlpha(k,n,e);
        
        % compute W2, W4
        [W2] = getW2(n,j,k,g);
        
        sumk1 = sumk1 + alp/2^k*func(k,j,n,i,e,eta)*W2;
        
    end
    
    if isinf(blp) == 0
        sumj1 = sumj1 + blp*sumk1;
    elseif sumk1 ~= 0
        disp('Inclination singularity');
    end
        
end

W1FW2 = deln*sumj1;

end

function [W1FW4] = getW1W4(func,n,Jn,J2,a,e,i,eta,g,mu,Re)

deln = deltalp(n,Jn,J2,a,eta,mu,Re);

sumj2 = 0;

% outer j loop
for j = 0:floor(n/2)
   
    blp = betalp(j,n,i);
    
    sumk2 = 0;
    
    % innner k loop
    for k1 = 0:(floor(n/2)-1)
        
        if mod(n,2) == 0
            % even harmonic
            k = 2*k1;
        else
            % odd harmonic
            k = 2*k1 + 1;
        end
        
        alp = getAlpha(k,n,e);
        
        % compute W2, W4
        [W4] = getW4(n,j,k,g);
        
        sumk2 = sumk2 + alp/2^k*func(k,j,n,i,e,eta)*W4;
    
    end

    if isinf(blp) == 0
        sumj2 = sumj2 + blp*sumk2;
    elseif sumk2 ~= 0
        disp('Inclination singularity');
    end    

end

W1FW4 = deln*sumj2;

end

function [W3FW2] = getW3W2(func,n,Jn,J2,a,e,i,eta,g,mu,Re)

deln = deltalp(n,Jn,J2,a,eta,mu,Re);

sumj1 = 0;

% outer j loop
for j = 0:floor(n/2)
   
    blp = betalp(j,n,i);
    
    sumk1 = 0;
    
    % innner k loop
    for k1 = 0:(floor(n/2)-1)
        
        if mod(n,2) == 0
            % even harmonic
            k = 2*k1;
        else
            % odd harmonic
            k = 2*k1 + 1;
        end
        
        if k < 1
            continue;
        end
        
        alp = getAlpha1(k,n,e);
        
        % compute W2, W4
        W2 = getW2(n,j,k,g);
        
        sumk1 = sumk1 + (n-k)/k*alp/2^k*func(k,j,n,i,e,eta)*W2;
    
    end

    if isinf(blp) == 0
        sumj1 = sumj1 + blp*sumk1;
    elseif sumk1 ~= 0
        disp('Inclination singularity');
    end    
        
end

W3FW2 = deln*sumj1;

end

function [W3FW2] = getW3bW2(func,n,Jn,J2,a,e,i,eta,g,mu,Re)

deln = deltalp(n,Jn,J2,a,eta,mu,Re);

sumj1 = 0;

% outer j loop
for j = 0:floor((n-2)/2)
   
    blp = betalp2(j,n,i);
    
    sumk1 = 0;
    
    % innner k loop
    for k1 = 0:(floor(n/2)-1)
        
        if mod(n,2) == 0
            % even harmonic
            k = 2*k1;
        else
            % odd harmonic
            k = 2*k1 + 1;
        end
        
        if k < 1
            continue;
        end
        
        alp = getAlpha1(k,n,e);
        
        % compute W2, W4
        W2 = getW2(n,j,k,g);
        
        sumk1 = sumk1 + (n-k)/k*alp/2^k*func(k,j,n,i,e,eta)*W2;
    
    end

    if isinf(blp) == 0
        sumj1 = sumj1 + blp*sumk1;
    elseif sumk1 ~= 0
        disp('Inclination singularity');
    end    
        
end

W3FW2 = deln*sumj1;

end

function [W3FW2] = getW3ebW2(func,n,Jn,J2,a,e,i,eta,g,mu,Re)

deln = deltalp(n,Jn,J2,a,eta,mu,Re);

sumj1 = 0;

% outer j loop
for j = 0:floor((n-2)/2)
   
    blp = betalp2(j,n,i);
    
    sumk1 = 0;
    
    % innner k loop
    for k1 = 1:(floor(n/2)-1)
        
        if mod(n,2) == 0
            % even harmonic
            k = 2*k1;
        else
            % odd harmonic
            k = 2*k1 + 1;
        end
        
        if k < 2
            continue;
        end
        
        alp = nchoosek(n-1,k-1)*e^(k-2);
        
        % compute W2, W4
        W2 = getW2(n,j,k,g);
        
        sumk1 = sumk1 + (n-k)/k*alp/2^k*func(k,j,n,i,e,eta)*W2;
    
    end

    if isinf(blp) == 0
        sumj1 = sumj1 + blp*sumk1;
    elseif sumk1 ~= 0
        disp('Inclination singularity');
    end    
        
end

W3FW2 = deln*sumj1;

end

function [W3FW2] = getW32W4(func,n,Jn,J2,a,e,i,eta,g,mu,Re)

deln = deltalp(n,Jn,J2,a,eta,mu,Re);

sumj1 = 0;

% outer j loop
for j = 0:floor((n-2)/2)
   
    blp = betalp2(j,n,i);
    
    sumk1 = 0;
    
    % innner k loop
    for k1 = 1:(floor(n/2)-1)
        
        if mod(n,2) == 0
            % even harmonic
            k = 2*k1;
        else
            % odd harmonic
            k = 2*k1 + 1;
        end
        
        if k < 2
            continue;
        end
        
        alp = nchoosek(n-1,k-1)*(k-1)*e^(k-2);
        
        % compute W2, W4
        W2 = getW4(n,j,k,g);
        
        sumk1 = sumk1 + (n-k)/k*alp/2^k*func(k,j,n,i,e,eta)*W2;
    
    end

    if isinf(blp) == 0
        sumj1 = sumj1 + blp*sumk1;
    elseif sumk1 ~= 0
        disp('Inclination singularity');
    end    
        
end

W3FW2 = deln*sumj1;

end

function [W3FW2] = getW3jstar(func,jstar,n,Jn,J2,a,e,i,eta,g,mu,Re)

if mod(n,2) == 0
    W3FW2 = 0;
else
    
    deln = deltalp(n,Jn,J2,a,eta,mu,Re);
    
    
    blp = betalp(jstar,n,i);
    
    sumk1 = 0;
    
    % innner k loop
    for k1 = 0:(floor(n/2)-1)
        
        if mod(n,2) == 0
            % even harmonic
            k = 2*k1;
        else
            % odd harmonic
            k = 2*k1 + 1;
        end
        
        if k < 1
            continue;
        end
        
        alp = getAlpha1(k,n,e);
        
        % get the customized function
        sumk1 = sumk1 + (n-k)/k*alp/2^k*func(k,jstar,n,i,e,eta);
        
    end
    
    W3FW2 = deln*blp*sumk1;
end
end

function [W3FW2] = getW3W2njp(func,n,Jn,J2,a,e,i,eta,g,mu,Re)

deln = deltalp(n,Jn,J2,a,eta,mu,Re);

sumj1 = 0;

% outer j loop
for j = 0:floor((n)/2)
   
    blp = betalp(j,n,i);
    
    sumk1 = 0;
    
    % innner k loop
    for k1 = 0:(floor(n/2)-1)
        
        if mod(n,2) == 0
            % even harmonic
            k = 2*k1;
        else
            % odd harmonic
            k = 2*k1 + 1;
        end
        
        if k < 1
            continue;
        end
        
        alp = getAlpha1(k,n,e);
        
        % compute W2, W4
        W2 = getW2njp(n,j,k,g);
        
        sumk1 = sumk1 + (n-k)/k*alp/2^k*func(k,j,n,i,e,eta)*W2;
    
    end

    if isinf(blp) == 0
        sumj1 = sumj1 + blp*sumk1;
    elseif sumk1 ~= 0
        disp('Inclination singularity');
    end    
        
end

W3FW2 = deln*sumj1;

end

function [W3FW2] = getW3bW2njp(func,n,Jn,J2,a,e,i,eta,g,mu,Re)

deln = deltalp(n,Jn,J2,a,eta,mu,Re);

sumj1 = 0;

% outer j loop
for j = 0:floor((n-2)/2)
   
    blp = betalp2(j,n,i);
    
    sumk1 = 0;
    
    % innner k loop
    for k1 = 0:(floor(n/2)-1)
        
        if mod(n,2) == 0
            % even harmonic
            k = 2*k1;
        else
            % odd harmonic
            k = 2*k1 + 1;
        end
        
        if k < 1
            continue;
        end
        
        alp = getAlpha1(k,n,e);
        
        % compute W2, W4
        W2 = getW2njp(n,j,k,g);
        
        sumk1 = sumk1 + (n-k)/k*alp/2^k*func(k,j,n,i,e,eta)*W2;
    
    end

    if isinf(blp) == 0
        sumj1 = sumj1 + blp*sumk1;
    elseif sumk1 ~= 0
        disp('Inclination singularity');
    end    
        
end

W3FW2 = deln*sumj1;

end

function [W3FW4] = getW3W4(func,n,Jn,J2,a,e,i,eta,g,mu,Re)

deln = deltalp(n,Jn,J2,a,eta,mu,Re);

sumj2 = 0;

% outer j loop
for j = 0:floor(n/2)
   
    blp = betalp(j,n,i);
    
    sumk2 = 0;
    
    % innner k loop
    for k1 = 0:(floor(n/2)-1)
        
        if mod(n,2) == 0
            % even harmonic
            k = 2*k1;
        else
            % odd harmonic
            k = 2*k1 + 1;
        end
        
        if k < 1
            continue;
        end
        
        alp = getAlpha1(k,n,e);
        
        % compute W2, W4
        W4 = getW4(n,j,k,g);
        
        sumk2 = sumk2 + (n-k)/k*alp/(2^k)*func(k,j,n,i,e,eta)*W4;
    
    end

    if isinf(blp) == 0
        sumj2 = sumj2 + blp*sumk2;
    elseif sumk2 ~= 0
        disp('Inclination singularity');
    end    
    
end

W3FW4 = deln*sumj2;

end

function [W3FW4] = getW3bW4(func,n,Jn,J2,a,e,i,eta,g,mu,Re)

deln = deltalp(n,Jn,J2,a,eta,mu,Re);

sumj2 = 0;

% outer j loop
for j = 0:floor((n-2)/2)
   
    blp = betalp2(j,n,i);
    
    sumk2 = 0;
    
    % innner k loop
    for k1 = 0:(floor(n/2)-1)
        
        if mod(n,2) == 0
            % even harmonic
            k = 2*k1;
        else
            % odd harmonic
            k = 2*k1 + 1;
        end
        
        if k < 1
            continue;
        end
        
        alp = getAlpha1(k,n,e);
        
        % compute W2, W4
        W4 = getW4(n,j,k,g);
        
        sumk2 = sumk2 + (n-k)/k*alp/(2^k)*func(k,j,n,i,e,eta)*W4;
    
    end

    if isinf(blp) == 0
        sumj2 = sumj2 + blp*sumk2;
    elseif sumk2 ~= 0
        disp('Inclination singularity');
    end    
    
end

W3FW4 = deln*sumj2;

end

function [W2] = getW2(n,j,k,g)

pend = floor((n - 1)/2) - j;

sumW2 = 0;

for p = 0:pend

    gam = getGamma(p,j,n);
    njp = n - 2*j - 2*p;
    
    if ((k - njp)/2) < 0 || ((k - njp)/2) > k
        continue;
    end
    
    bc = nchoosek(k, (k - njp)/2 );
    
    sgt = sin(njp*g);
    cgt = cos(njp*g);

    if mod(n,2) == 0
  
        % even harmonic
        sumW2 = sumW2 + gam*bc*(-sgt)/njp;
    else
        % odd harmonic
        sumW2 = sumW2 + gam*bc*(cgt)/njp;
    end
    
end

W2 = sumW2;

end

function [W2] = getW2njp(n,j,k,g)

pend = floor((n - 1)/2) - j;

sumW2 = 0;

for p = 0:pend

    gam = getGamma(p,j,n);
    njp = n - 2*j - 2*p;
    
    if ((k - njp)/2) < 0 || ((k - njp)/2) > k
        continue;
    end
    
    bc = nchoosek(k, (k - njp)/2 );
    
    sgt = sin(njp*g);
    cgt = cos(njp*g);

    if mod(n,2) == 0
  
        % even harmonic
        sumW2 = sumW2 + gam*bc*(-sgt)/njp*njp^2;
    else
        % odd harmonic
        sumW2 = sumW2 + gam*bc*(cgt)/njp*njp^2;
    end
    
end

W2 = sumW2;

end

function [W4] = getW4(n,j,k,g)

pend = floor((n - 1)/2) - j;

sumW4 = 0;

for p = 0:pend

    gam = getGamma(p,j,n);
    njp = n - 2*j - 2*p;
    
    if ((k - njp)/2) < 0 || ((k - njp)/2) > k
        continue;
    end
    
    bc = nchoosek(k, (k - njp)/2 );
    
    sgt = sin(njp*g);
    cgt = cos(njp*g);

    
    if mod(n,2) == 0
  
        % even harmonic
        sumW4 = sumW4 + gam*bc*(-cgt);
    else
        % odd harmonic
        sumW4 = sumW4 + gam*bc*(-sgt);
    end
    
end

W4 = sumW4;

end

function gam = getGamma(p,j,n)

s = 1;
if mod(n,2) ~= 0 
    % odd harmonic
    s = -1;
end

gam = ((-1)^(floor(n/2) - j - (s*p)))*2*nchoosek(n - 2*j, p);

end

function dellp = deltalp(n,Jn,J2,a,eta,mu,Re)
dellp = 4/3*Jn/J2^2*sqrt(mu)*Re^(n-2)*a^(7/2-n-1)*eta^(4-2*n+1)/2^n;
end

function blp = betalp(j,n,i)

blp = (-1)^j*factorial(2*n-2*j)*(sin(i))^(n-2*j-1)/...
        ((1-5*(cos(i))^2)*factorial(j)*factorial(n-j)*factorial(n-2*j)*2^(n-2*j));

end

function blp = betalp2(j,n,i)

blp = (-1)^j*factorial(2*n-2*j)*(sin(i))^(n-2*j-2)/...
        ((1-5*(cos(i))^2)*factorial(j)*factorial(n-j)*factorial(n-2*j)*2^(n-2*j));

end

function blp = betalpjf(j,n,i,jf)

blp = (-1)^j*factorial(2*n-2*j)*(sin(i))^(n-2*j-jf)/...
        ((1-5*(cos(i))^2)*factorial(j)*factorial(n-j)*factorial(n-2*j)*2^(n-2*j));

end

function alp = getAlpha(k,n,e)

alp = e^k*nchoosek(n-1,k);

end

function alp = getAlpha1(k,n,e)

alp = e^(k-1)*nchoosek(n-1,k-1);

end

