%% EqSPZonal
%     Desc: Computes 2nd order Short-period effects for Jn (n >= 3) for Equinoctial elements.
%     Equinoctial elements: [a l+g+h tan(i/2)*cos(h) tan(i/2)*sin(h) e*cos(g+h) e*sin(g+h)];
% Reference: 
% B. Mahajan, S. R. Vadali and K. T. Alfriend, �Analytic Solution For Satellite 
% Relative Motion: The Complete Zonal Gravitational Problem,� 26th AAS/AIAA 
% Space Flight Mechanics Meeting, Napa, CA, February 2016.
% Author: Bharat Mahajan (https://github.com/princemahajan)

function [EqSP, Dsp] = EqSPZonal(Xm,degree,mu,Re,Jcoeff,InverseOn, JacobianOn, tol) 

% The following value will be returned by beta, betaSP functions in case
% the divide by 0 singularity occurs in these functions.
ERR_VAL = -999999999;

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

% intermediate quantities
psi = Mean2TrueLong(Lambda, q1, q2, tol);
lop = mod(atan2(q2,q1),2*pi); 
f = psi - lop;
Psi = psi;
Delta = atan((tan(Psi)-tan(Lambda))/(1+tan(Psi)*tan(Lambda)));

dfdl = (1 + e*cos(f))^2/eta^3;

dfde = (2+e*cos(f))*sin(f)/eta^2;

d2fdfdlbye = -2*(1+e*cos(f))*sin(f)/eta*3;

d2fdfdl = e*d2fdfdlbye;

d2fdfdldfdlbye = 1/2/eta^6*((3*eta^2-7)*sin(f) + e*(eta^2-7)*sin(2*f) + (3*eta^2-3)*sin(3*f) ...
                    + e/2*(eta^2-1)*sin(4*f));

d2fdl2bye = d2fdfdldfdlbye;

d2fdl2 = e*d2fdl2bye;

dfdlm1bye = 2*cos(f)/eta^3 + 1/2*e*cos(2*f)/eta^3 + 1/2*e*(2*eta^2 + 3*eta + 3)/(1+eta)/eta^3;

d2fdlde = 1/2/eta^5*((9-5*eta^2)*cos(f) - e*(eta^2-7)*cos(2*f) + (3-3*eta^2)*cos(3*f) ...
            - e/2*(eta^2-1)*cos(4*f) - e/2*(eta^2 - 9));

fLG = -sqrt((1-eta)/(1+eta))*sin(f)*(2+e*cos(f))/G;

fLGdfbye = (-2*cos(f)-e*cos(2*f))/L/eta/(1+eta);

fLGdf = fLGdfbye*e;

dfLGde = 1/2/sqrt(mu*a)/eta^3*((3*eta^2-4*eta-3)*sin(f)/(1+eta) - e*(eta+5)*sin(2*f)/(1+eta) ...
            + 3*(eta-1)*sin(3*f) + 1/2*(eta-1)*e*sin(4*f));

fLGbye = 1/L/(1+eta)/eta*(-2*sin(f) - 1/2*e*sin(2*f));        

dfdlfLGbye = dfdl*fLGbye;

d2fdfde = -e * sin(f) ^ 2 / eta ^ 2 + (0.2e1 + e * cos(f)) * cos(f) / eta ^ 2;

dfdee = 1/eta^4*(4*e*sin(f) + 1/2*(-eta^2 + 2)*sin(2*f));

d2fde2 = d2fdfde*dfde + dfdee;
        
% J2 Equinoctial First Order Short-Period Effects (Using Maple)

% J2 Equinoctial Second Order Short-Period Effects (Using Maple)


% Second order Short-Period effects due to higher degree zonals
EqSP = zeros(6,1);
Dsp = zeros(6);

for ctr = 1:length(narr)
    
    n = narr(ctr);
    Jn = Jcoeff(n);

    % Short-period Effects on Equinoctial elements

    % intermediate terms
    [K2n,K2nk1bye,K2nk2bye] = getK2nF(@dummy, n,Jn,J2,a,e,eta,i,g,G,f,mu,Re,false);
    
    [K2nLGH,K2nLGHk1bye] = getK2nF(@KLGHfunc, n,Jn,J2,a,e,eta,i,g,G,f,mu,Re,false);
    
    [K2nn2j,~,~] = getK2nF(@n2jfunc, n,Jn,J2,a,e,eta,i,g,G,f,mu,Re,false);
    K2nn2jbysi = getK2nF(@n2jfunc,  n,Jn,J2,a,e,eta,i,g,G,f,mu,Re,true);
    
    [K2ngbysi,K2ngk1byebysi,~] = getK2ngF(@dummy,n,Jn,J2,a,e,eta,i,g,G,f,mu,Re,true);
    K2ng = K2ngbysi*si;
    K2ngk1bye = K2ngk1byebysi*si;
        
    K2neGH1 = getK2nF(@KeGHfunc1, n,Jn,J2,a,e,eta,i,g,G,f,mu,Re,false);
    [~,K2neGH2,~] = getK2nF(@KeGHfunc2, n,Jn,J2,a,e,eta,i,g,G,f,mu,Re,false);
    K2neGH = K2neGH1 + K2neGH2;

    
    W1 = sqrt(a^3/mu)*Delta*K2n;    

    
    [W2W3f, W2W3fk1bye, W2W3fk2bye] = getW2W3f(@dummy,n,Jn,J2,a,e,eta,i,f,l,G,g,mu,Re);
        
    W2n2jW3bysi = getW2W3F(@n2jfunc,n,Jn,J2,a,e,eta,i,f,l,G,g,mu,Re,true);
    W2n2jW3 = W2n2jW3bysi*si;

    W2W3gbysi = getW2W3g(@dummy,n,Jn,J2,a,e,eta,i,f,l,G,g,mu,Re,true);
    W2W3g = W2W3gbysi*si;

    W2eW3 = getW2eW3F(@dummy,n,Jn,J2,a,e,eta,i,f,l,G,g,mu,Re,false);    

    W2iW3 = getW2iW3(@dummy,n,Jn,J2,a,e,eta,i,f,l,G,g,mu,Re,false);    

    [W1lge, W1lgel] = getW1lge(n,Jn,J2,a,e,eta,i,f,l,g,G,Delta,mu,Re,K2n, K2ngk1bye);
    
    W23lge = getW23lge(@dummy,n,Jn,J2,a,e,eta,i,f,l,g,L,G,mu,Re,false);
    
    
    % SMA
    ansp2 = -2*sqrt(a/mu)*(sqrt(a^3/mu)*(dfdl-1)*K2n + W2W3f*dfdl);

    % MAOL
    [W2Lt3W3,~,~] = getW2W3F(@Lt3func,n,Jn,J2,a,e,eta,i,f,l,G,g,mu,Re,false);
    
    Lambdasp2 = 3/L*W1 + sqrt(a^3/mu)*(Delta*K2nLGH + K2n*fLG) - W2Lt3W3 + W2W3f*fLG;
    
    % Inclination terms p1 p2

    if ci == 0
        p1sp2 = 0; 
        p2sp2 = 0;
    else
        t1p12 = ci/G/(1+ci);
        t2p12o = sqrt(a^3/mu)*Delta*(-K2nn2jbysi) - W2n2jW3bysi;
        t3p12o = sqrt(a^3/mu)*Delta*K2ngbysi + W2W3gbysi;
        
        p1sp2 = t1p12*(-sh*t2p12o - ch*t3p12o);
        p2sp2 = t1p12*(ch*t2p12o - sh*t3p12o);
    end
    
    % Eccentricity terms q1 q2
    
    t1q12o = sqrt(a^3/mu)*(Delta*K2neGH - dfde*eta/L*K2n);
    t2q12o = -eta/L*W2eW3 - e/G/(1+ci)*W2iW3 - eta/L*dfde*W2W3f;
    t3q12o = W1lge + W23lge;
    
    q1sp2 = -sin(g+h)*(t1q12o + t2q12o) - cos(g+h)*t3q12o;
    q2sp2 = cos(g+h)*(t1q12o + t2q12o) - sin(g+h)*t3q12o;
    
    JnSP = [ansp2, Lambdasp2, p1sp2,p2sp2,q1sp2,q2sp2]';
    EqSP = EqSP + JnSP;

    if JacobianOn ~= true
        continue;
    end
    
    %  Generate D matrix

    % intermediate terms

    [K2ngbysi,K2ngk1byebysi,~] = getK2ngF(@dummy,n,Jn,J2,a,e,eta,i,g,G,f,mu,Re,true);
    K2ngbye = K2ngk1byebysi*si;
    [~,K2nk1byek,~] = getK2nF(@kfunc, n,Jn,J2,a,e,eta,i,g,G,f,mu,Re,false);
    K2ne = (2*n-1)*e/eta^2*K2n + K2nk1byek;
    [~,K2nk1byekn2j,~] = getK2nF(@kn2jfunc, n,Jn,J2,a,e,eta,i,g,G,f,mu,Re,false);
    K2nen2j = (2*n-1)*e/eta^2*K2nn2j + K2nk1byekn2j;
    [~,K2ngk1byek,~] = getK2ngF(@kfunc, n,Jn,J2,a,e,eta,i,g,G,f,mu,Re,false);
    K2nge = (2*n-1)*e/eta^2*K2ng + K2ngk1byek;    
    dK2ndi = ci*K2nn2jbysi;
    
    K2nLGHn2jbysi = getK2nF(@KLGHn2jfunc, n,Jn,J2,a,e,eta,i,g,G,f,mu,Re,true);
    K2nLGHdi2 = getK2nF(@KLGHdi2func, n,Jn,J2,a,e,eta,i,g,G,f,mu,Re,false);
    dK2nLGHdi = ci*K2nLGHn2jbysi + K2nLGHdi2;

    K2nn2jbysi = getK2nF(@n2jfunc, n,Jn,J2,a,e,eta,i,g,G,f,mu,Re,true);    
    [K2nLGHgbysi,K2nLGHgk1bysibye,~] = getK2ngF(@KLGHfunc,n,Jn,J2,a,e,eta,i,g,G,f,mu,Re,true);
    
    [~,K2nLGHk1byek] = getK2nF(@KLGHkfunc, n,Jn,J2,a,e,eta,i,g,G,f,mu,Re,false);
    K2nLGHe = (2*n-1)*e/eta^2*K2nLGH + e/eta^2*K2nLGH + K2nLGHk1byek;
    K2nn2jn2jbysi = getK2nF(@n2jn2jfunc, n,Jn,J2,a,e,eta,i,g,G,f,mu,Re,true);    
    [K2ngn2jbysi,K2ngn2jk1byebysi,~] = getK2ngF(@n2j1func,n,Jn,J2,a,e,eta,i,g,G,f,mu,Re,true);
    [K2ng2bysi,K2ng2k1byebysi,~] = getK2ng2F(@dummy, n,Jn,J2,a,e,eta,i,g,G,f,mu,Re,true);
    K2ngn2jbye = K2ngn2jk1byebysi*si;
        
    K2neGH1ibysi = getK2nF(@KeGHfunc1i, n,Jn,J2,a,e,eta,i,g,G,f,mu,Re,true);
    [~,K2neGH2ibysi,~] = getK2nF(@KeGHfunc2i, n,Jn,J2,a,e,eta,i,g,G,f,mu,Re,true);
    K2neGH3i = getK2nF(@KeGHfunc3i, n,Jn,J2,a,e,eta,i,g,G,f,mu,Re,false);
    K2neGHi = (K2neGH1ibysi + K2neGH2ibysi)*si + K2neGH3i;
    
    K2neGHg1bysi = getK2ngF(@KeGHfunc1, n,Jn,J2,a,e,eta,i,g,G,f,mu,Re,true);
    [~,K2neGHg2k1byebysi] = getK2ngF(@KeGHgfunc2, n,Jn,J2,a,e,eta,i,g,G,f,mu,Re,true);
    K2neGHgbysi = K2neGHg1bysi + K2neGHg2k1byebysi;  
    
    K2nSPeGHe1 = getK2nF(@KeGHefunc1, n,Jn,J2,a,e,eta,i,g,G,f,mu,Re,false);
    [~,K2neGHek1bye,~] = getK2nF(@KeGHefunc2, n,Jn,J2,a,e,eta,i,g,G,f,mu,Re,false);
    [~,~,K2neGHek2bye] = getK2nF(@KeGHefunc3, n,Jn,J2,a,e,eta,i,g,G,f,mu,Re,false);    
    K2neGHe = (2*n*e/eta^2)*K2neGH + K2nSPeGHe1 + K2neGHek1bye + K2neGHek2bye;        
    
    W1lbye = sqrt(a^3/mu)*K2n*dfdlm1bye;
    W1l = W1lbye*e;
    W1n2jbysi = sqrt(a^3/mu)*Delta*K2nn2jbysi;
    W1gbysi = sqrt(a^3/mu)*Delta*K2ngbysi;
    W1e = sqrt(a^3/mu)*(dfde*K2n + Delta*K2ne);
    W1gbye = sqrt(a^3/mu)*(f-l)*K2ngk1bye;
    W1lgei = getW1lgei(n,Jn,J2,a,e,eta,i,f,l,Delta, g,G,mu,Re,dK2ndi);
    [W1lgegbysi] = getW1lgeg(n,Jn,J2,a,e,eta,i,f,l,Delta,g,G,mu,Re,K2ngbysi);
    [W1lgee] = getdW1lgede(n,Jn,J2,a,e,eta,i,f,l,Delta,g,G,mu,Re,K2n,K2ne,dfde);

    W2iW3fbysi = getW2iW3f(@dummy,n,Jn,J2,a,e,eta,i,f,l,G,g,mu,Re,true);    
    W2iW3f = W2iW3fbysi*si;
    [W2W3fgbysi,W2W3fgk1byebysi,W2W3fgk2byebysi] = getW2W3fg(@dummy,n,Jn,J2,a,e,eta,i,f,l,G,g,mu,Re,true);
    W2W3fg = W2W3fgbysi*si;
    W2W3fgk1bye = W2W3fgk1byebysi*si;
    W2W3fgk2bye = W2W3fgk2byebysi*si;
    W2eW3f = getW2eW3fF(@dummy,n,Jn,J2,a,e,eta,i,f,l,G,g,mu,Re,false);
    [W2W3f2,W2W3f2k1bye,W2W3f2k2bye] = getW2W3f2(@dummy,n,Jn,J2,a,e,eta,i,f,l,G,g,mu,Re,false);
    [W2W3fn2jbysi,W2W3fn2jk1byebysi,~] = getW2W3f(@n2jfunc,n,Jn,J2,a,e,eta,i,f,l,G,g,mu,Re);
    W2W3fn2j = W2W3fn2jbysi*si;
    W2W3fn2jk1bye = W2W3fn2jk1byebysi*si;
    W2n2jn2jW3bysi = getW2W3F(@n2jn2jfunc,n,Jn,J2,a,e,eta,i,f,l,G,g,mu,Re,true);
    [W2n2jW3gbysi,W2n2jW3gk1byebysi,~] = getW2W3g(@n2jfunc,n,Jn,J2,a,e,eta,i,f,l,G,g,mu,Re,true);
    [W2W3g2bysi,W2W3g2k1byebysi,W2W3g2k2byebysi] = getW2W3g2(@dummy,n,Jn,J2,a,e,eta,i,f,l,G,g,mu,Re,true);
    W2W3g2k2bye = W2W3g2k2byebysi*si;
    W2en2jW3 = getW2eW3F(@n2jfunc,n,Jn,J2,a,e,eta,i,f,l,G,g,mu,Re,false);
    
    W23lgel = getW23lgel(n,Jn,J2,a,e,eta,i,f,l,g,L,G,mu,Re,dfdl,d2fdfdl,false);
    W2eiW3 = ci*getW2eW3F(@n2jfunc,n,Jn,J2,a,e,eta,i,f,l,G,g,mu,Re,true);
    W2i2W3 = W2n2jn2jW3bysi*ci^2;
        
    W23lgei = ci*getW23lge(@n2jfunc,n,Jn,J2,a,e,eta,i,f,l,g,L,G,mu,Re,true);
    W2iW3gbysi = getW2iW3g(@dummy,n,Jn,J2,a,e,eta,i,f,l,G,g,mu,Re,true);
    W2iW3g = W2iW3gbysi*si;
    W2eW3gbysi = getW2eW3gF(@dummy,n,Jn,J2,a,e,eta,i,f,l,G,g,mu,Re,true);
    W2eW3g = W2eW3gbysi*si;
    W23lgegbysi = getW23lgeg(@dummy,n,Jn,J2,a,e,eta,i,f,l,g,L,G,dfdl,mu,Re);
    W2e2W3 = getW2e2W3F(@dummy,n,Jn,J2,a,e,eta,i,f,l,G,g,mu,Re);
    [W2W3,~,~] = getW2W3F(@dummy,n,Jn,J2,a,e,eta,i,f,l,G,g,mu,Re,false);
    [~,~,W2W3k2byek] = getW2W3F(@kterm,n,Jn,J2,a,e,eta,i,f,l,G,g,mu,Re,false);
    W23lgee = getdW23lgede(n,Jn,J2,a,e,eta,i,f,l,g,L,G,W23lge,mu,Re,dfdl,dfde,d2fdlde);
        
    
%     K2nSPgbye = getK2nSPgbyeF(@dummy,n,Jn,J2,a,e,eta,i,g,mu,Re);
%     K2nSPge = (2*n-1)*e/eta^2*K2nSPg + getK2nSPgF(@kfunc,n,Jn,J2,a,e,eta,i,g,mu,Re);
%     [K2nSPg2,K2nSPg2bye] = getK2nSPg2F(@dummy,n,Jn,J2,a,e,eta,i,g,mu,Re);
%     W2eW3g = getW2eW3g(n,Jn,J2,a,e,eta,i,f,l,G,g,mu,Re);
%     [W2n2jW3g,W2n2jW3gk1bye] = getW2W3g(@n2jfunc,n,Jn,J2,a,e,eta,i,f,l,G,g,mu,Re);
%     W2W3g2m1 = getW2W3g2m1F(@dummy,n,Jn,J2,a,e,eta,i,f,l,G,g,mu,Re);
%     W2iW3gSI = getW2iW3gSI(@dummy,n,Jn,J2,a,e,eta,i,f,l,G,g,mu,Re);
    
    % SMA Row
    
    % da
    daspda = -(n-1)/a*ansp2;
    
    % dL
    daspdL = -2*sqrt(a/mu)*(sqrt(a^3/mu)*K2n*d2fdl2 + W2W3f*d2fdl2 + W2W3f2*(dfdl)^2);
    
    % dp1 dp2
    dW2diW3f = W2iW3fbysi;
    daspdi = -2*sqrt(a/mu)*(sqrt(a^3/mu)*(dfdl-1)*dK2ndi + dW2diW3f*dfdl);
    
    daspdgbysi = -2*sqrt(a/mu)*(sqrt(a^3/mu)*(dfdl-1)*K2ngbysi + W2W3fgbysi*dfdl);
    
    daspdp1 = daspdi*(1+ci)*ch + daspdgbysi*(1+ci)*sh;
    
    daspdp2 = daspdi*(1+ci)*sh - daspdgbysi*(1+ci)*ch;
    
    % dq1 dq2
    daspde = -2*sqrt(a/mu)*(sqrt(a^3/mu)*(dfdl-1)*K2ne + sqrt(a^3/mu)*K2n*d2fdlde ...
                + W2eW3f*dfdl + W2W3f2*dfde*dfdl + W2W3f*d2fdlde);
            
    daspdlq12 = -2*sqrt(a/mu)*(sqrt(a^3/mu)*d2fdl2bye*K2n + W2W3f*d2fdl2bye + W2W3f2k1bye*dfdl^2);
    
    daspdgq12 = -2*sqrt(a/mu)*(sqrt(a^3/mu)*dfdlm1bye*K2ng + W2W3fgk1bye*dfdl);
    
    W2k0W3f2fg = -2*sqrt(a/mu)*dfdl*dfdlm1bye*getW2k0F(@k0dadldadgterm,n,Jn,J2,a,e,eta,i,f,l,G,g,mu,Re);
    
    daspdq1 = daspde*cos(g+h) + daspdlq12*sin(g+h) - daspdgq12*sin(g+h) + W2k0W3f2fg*sin(g+h);
    
    daspdq2 = daspde*sin(g+h) - daspdlq12*cos(g+h) + daspdgq12*cos(g+h) - W2k0W3f2fg*cos(g+h);
    
    % Mean Longitude Row
    [W2Lt3W3f,W2Lt3W3fk1bye] = getW2W3f(@Lt3func,n,Jn,J2,a,e,eta,i,f,l,G,g,mu,Re);
    W2Lt3W3n2jbysi = getW2W3F(@Lt3n2jfunc,n,Jn,J2,a,e,eta,i,f,l,G,g,mu,Re,true);
    W2Lt3W3i2 = getW2W3F(@Lt3i2func,n,Jn,J2,a,e,eta,i,f,l,G,g,mu,Re,false);
    [W2Lt3W3gbysi,W2Lt3W3gk1byebysi] = getW2W3g(@Lt3func,n,Jn,J2,a,e,eta,i,f,l,G,g,mu,Re,true);
    W2eLt3W3 = getW2eW3F(@Lt3func,n,Jn,J2,a,e,eta,i,f,l,G,g,mu,Re,false);
    W2Lt3W3e2 = getW2W3F(@Lt3e2func,n,Jn,J2,a,e,eta,i,f,l,G,g,mu,Re,false);

    % da
    dLspda = 3/L/2/a*W1*(-2*n) + Delta*K2nLGH/2*sqrt(a/mu)*(-2*n) + K2n*fLG*sqrt(a/mu)*(-n) ...
                    + 2*n/2/a*W2Lt3W3 - (2*n-1)/2/a*W2W3f*fLG - 1/2/a*W2W3f*fLG;
    
    % dLambda
    dLspdL = 3/L*W1l + sqrt(a^3/mu)*(dfdlm1bye*e*K2nLGH + K2n*fLGdf*dfdl) - W2Lt3W3f*dfdl  ...
                + W2W3f2*dfdl*fLG + W2W3f*dfdl*fLGdf;
    
    % dp1 dp2
    dLspdi = 3/L*W1n2jbysi*ci + sqrt(a^3/mu)*(Delta*dK2nLGHdi + K2nn2jbysi*ci*fLG) ...
                - W2Lt3W3n2jbysi*ci - W2Lt3W3i2 + W2W3fn2jbysi*ci*fLG;
   
    dLspdgbysi = 3/L*W1gbysi + sqrt(a^3/mu)*(Delta*K2nLGHgbysi + K2ngbysi*fLG) ...
                    - W2Lt3W3gbysi + W2W3fgbysi*fLG;
            
    dLspdp1 = dLspdi*ch*(1+ci) + dLspdgbysi*sh*(1+ci);
    
    dLspdp2 = dLspdi*sh*(1+ci) - dLspdgbysi*ch*(1+ci);
    
    % dq1 dq2
    
    dLspde = 3/L*W1e + sqrt(a^3/mu)*(Delta*K2nLGHe + dfde*K2nLGHe + K2ne*fLG + K2n*dfLGde) ...
            - W2eLt3W3 - W2Lt3W3*e/eta^2 - W2Lt3W3e2 + W2eW3f*fLG + W2W3f2*dfde*fLG + W2W3f*dfLGde;

    dLspdl = 3/L*W1lbye + sqrt(a^3/mu)*(dfdlm1bye*K2nLGH + K2n*fLGdfbye*dfdl) - W2Lt3W3fk1bye*dfdl  ...
                + W2W3f2*dfdlfLGbye + W2W3f*dfdl*fLGdfbye;
    
    dLspdg = 3/L*W1gbye + sqrt(a^3/mu)*(Delta*K2nLGHgk1bysibye*si + K2ngbye*fLG)...
                - W2Lt3W3gk1byebysi*si + W2W3fg*fLGbye;

    W2Lt3W3fgk0 = -getW2k0F(@k0W2Lt3W3fgterm,n,Jn,J2,a,e,eta,i,f,l,G,g,mu,Re)*dfdlm1bye;            
            
    dLspdq1 = dLspde*cos(g+h) + dLspdl*sin(g+h) - dLspdg*sin(g+h) + W2Lt3W3fgk0*sin(g+h);
            
    dLspdq2 = dLspde*sin(g+h) - dLspdl*cos(g+h) + dLspdg*cos(g+h) - W2Lt3W3fgk0*cos(g+h);
    
    % Incination row p1,p2
    
    % da
    dp1spda = -n/a*(p1sp2);
    dp2spda = -n/a*(p2sp2);

    % dL
    t1 = ci/G/(1+ci);
    t2 = sqrt(a^3/mu)*dfdlm1bye*e*(-K2nn2j) - W2W3fn2jbysi*si*dfdl;
    t3 = sqrt(a^3/mu)*dfdlm1bye*e*K2ng + W2W3fg*dfdl;
    dp1spdL = t1*(-sh*t2 - ch*t3);
    dp2spdL = t1*(ch*t2 - sh*t3);
    
    % dp1 dp2
    dp1spdi = -si/G/(1+ci)^2*(-sh*t2p12o - ch*t3p12o) ...
                + t1p12*(-sh*(sqrt(a^3/mu)*Delta*(-K2nn2jn2jbysi*ci) - W2n2jn2jW3bysi*ci) ...
                -ch*(sqrt(a^3/mu)*Delta*K2ngn2jbysi*ci + W2n2jW3gbysi*ci) );
    
    dp1spdg = t1p12*(-sh*(sqrt(a^3/mu)*Delta*(-K2ngn2jbysi) - W2n2jW3gbysi) ...
                -ch*(sqrt(a^3/mu)*Delta*K2ng2bysi + W2W3g2bysi));
    
    dp1spdh = t1p12*(-ch*(sqrt(a^3/mu)*Delta*(-K2nn2jbysi) - W2n2jW3bysi) ...
                + sh*(sqrt(a^3/mu)*Delta*K2ngbysi) + W2W3gbysi);
            
    dp1spdp1 = dp1spdi*ch*(1+ci) + dp1spdg*sh*(1+ci) - dp1spdh*sh*(1+ci);
    
    dp1spdp2 = dp1spdi*sh*(1+ci) - dp1spdg*ch*(1+ci) + dp1spdh*ch*(1+ci);
    
    dp2spdi = -si/G/(1+ci)^2*(ch*t2p12o - sh*t3p12o) ...
                + t1p12*(ch*(sqrt(a^3/mu)*Delta*(-K2nn2jn2jbysi*ci) - W2n2jn2jW3bysi*ci) ...
                -sh*(sqrt(a^3/mu)*Delta*K2ngn2jbysi*ci + W2n2jW3gbysi*ci) );
    
    dp2spdg = t1p12*(ch*(sqrt(a^3/mu)*Delta*(-K2ngn2jbysi) - W2n2jW3gbysi) ...
                -sh*(sqrt(a^3/mu)*Delta*K2ng2bysi + W2W3g2bysi));
    
    dp2spdh = t1p12*(-sh*(sqrt(a^3/mu)*Delta*(-K2nn2jbysi) - W2n2jW3bysi) ...
                - ch*(sqrt(a^3/mu)*Delta*K2ngbysi) + W2W3gbysi);
            
    dp2spdp1 = dp2spdi*ch*(1+ci) + dp2spdg*sh*(1+ci) - dp2spdh*sh*(1+ci);
    
    dp2spdp2 = dp2spdi*sh*(1+ci) - dp2spdg*ch*(1+ci) + dp2spdh*ch*(1+ci);
    
    % dq1 dq2
    dp1spde = t1p12*e/eta^2*(-sh*t2p12o - ch*t3p12o) ...
                + t1p12*(-sh*(sqrt(a^3/mu)*dfde*(-K2nn2j) + sqrt(a^3/mu)*Delta*(-K2nen2j) ...
                - W2en2jW3 - W2W3fn2j*dfde) -ch*(sqrt(a^3/mu)*dfde*K2ng ...
                + sqrt(a^3/mu)*Delta*K2nge + W2eW3g + W2W3fg*dfde));
    
    dp1spdg = t1p12*(-sh*(sqrt(a^3/mu)*Delta*(-K2ngn2jbye) - W2n2jW3gk1byebysi*si) ...
                -ch*(sqrt(a^3/mu)*Delta*K2ng2k1byebysi*si + W2W3g2k1byebysi*si));
    
    dp1spdl = t1p12*(-sh*(sqrt(a^3/mu)*(dfdlm1bye)*(-K2nn2j) - W2W3fn2jk1bye*dfdl) ...
                -ch*(sqrt(a^3/mu)*(dfdlm1bye)*K2ng + W2W3fgk1byebysi*si*dfdl));
            
    dp1spdlgk0term = t1p12*getW2k0F(@dp1spdlgk0func,n,Jn,J2,a,e,eta,i,f,l,G,g,mu,Re);
    dp1spdlg2k0term = t1p12*getW2k0F(@dp1spdlg2k0func,n,Jn,J2,a,e,eta,i,f,l,G,g,mu,Re);
            
    dp1spdq1 = dp1spde*cos(g+h) - dp1spdg*sin(g+h) + dp1spdl*sin(g+h) ...
                + dp1spdlgk0term*sh*dfdlm1bye*sin(g+h) + dp1spdlg2k0term*ch*dfdlm1bye*sin(g+h);
            
    dp1spdq2 = dp1spde*sin(g+h) + dp1spdg*cos(g+h) - dp1spdl*cos(g+h) ...
                - dp1spdlgk0term*sh*dfdlm1bye*cos(g+h) - dp1spdlg2k0term*ch*dfdlm1bye*cos(g+h);            

    dp2spde = t1p12*e/eta^2*(ch*t2p12o - sh*t3p12o) ...
                + t1p12*(ch*(sqrt(a^3/mu)*dfde*(-K2nn2j) + sqrt(a^3/mu)*Delta*(-K2nen2j) ...
                - W2en2jW3 - W2W3fn2j*dfde) -sh*(sqrt(a^3/mu)*dfde*K2ng ...
                + sqrt(a^3/mu)*Delta*K2nge + W2eW3g + W2W3fg*dfde));
    
    dp2spdg = t1p12*(ch*(sqrt(a^3/mu)*Delta*(-K2ngn2jbye) - W2n2jW3gk1byebysi*si) ...
                -sh*(sqrt(a^3/mu)*Delta*K2ng2k1byebysi*si + W2W3g2k1byebysi*si));
    
    dp2spdl = t1p12*(ch*(sqrt(a^3/mu)*(dfdlm1bye)*(-K2nn2j) - W2W3fn2jk1bye*dfdl) ...
                - sh*(sqrt(a^3/mu)*(dfdlm1bye)*K2ng + W2W3fgk1byebysi*si*dfdl));
            
    dp2spdlgk0term = t1*getW2k0F(@dp1spdlgk0func,n,Jn,J2,a,e,eta,i,f,l,G,g,mu,Re);
    dp2spdlg2k0term = t1*getW2k0F(@dp1spdlg2k0func,n,Jn,J2,a,e,eta,i,f,l,G,g,mu,Re);
            
    dp2spdq1 = dp2spde*cos(g+h) - dp2spdg*sin(g+h) + dp2spdl*sin(g+h) ...
                - dp2spdlgk0term*ch*dfdlm1bye*sin(g+h) + dp2spdlg2k0term*sh*dfdlm1bye*sin(g+h);
            
    dp2spdq2 = dp2spde*sin(g+h) + dp2spdg*cos(g+h) - dp2spdl*cos(g+h) ...
                + dp2spdlgk0term*sh*dfdlm1bye*cos(g+h) - dp2spdlg2k0term*ch*dfdlm1bye*cos(g+h);            
            
    % Eccentricity terms row q1,q2
    
    % da
    dq1spda = -n/a*q1sp2;
    dq2spda = -n/a*q2sp2;

    % dL
    t1dl = sqrt(a^3/mu)*(dfdlm1bye*e*(K2neGH) - d2fdlde*eta/L*K2n);
    t2dl = -eta/L*W2eW3f*dfdl - e/G/(1+ci)*W2iW3fbysi*si*dfdl - eta/L*d2fdlde*W2W3f - eta/L*dfde*W2W3f2*dfdl;
    t3dl = W1lgel + W23lgel;
    
    dq1spdL = -sin(g+h)*(t1dl + t2dl) - cos(g+h)*t3dl;
    dq2spdL = cos(g+h)*(t1dl + t2dl) - sin(g+h)*t3dl;

    % dp1 dp2
    t1isi = sqrt(a^3/mu)*(Delta*(K2neGHi) - dfde*eta/L*K2nn2jbysi*ci);
    t2isi = -eta/L*W2eiW3 - e/G/(1+ci)*W2i2W3 + e/G*W2iW3*si/(1+ci)^2 - eta/L*dfde*W2iW3fbysi;
    t3isi = W1lgei + W23lgei;

    dq1spdi = - sin(g+h)*(t1isi + t2isi) - cos(g+h)*t3isi;
    
    t1g = sqrt(a^3/mu)*(Delta*(K2neGHgbysi) - dfde*eta/L*K2ngbysi);
    t2g = -eta/L*W2eW3gbysi - e/G/(1+ci)*W2iW3gbysi - eta/L*dfde*W2W3fgbysi;
    t3g = W1lgegbysi + W23lgegbysi;
    dq1spdgbysi = (-cos(g+h)*(t1q12o + t2q12o) + sin(g+h)*t3q12o)*0 + (-sin(g+h)*(t1g + t2g) - cos(g+h)*t3g);

    dq1spdhbysi = (-cos(g+h)*(t1q12o + t2q12o) + sin(g+h)*t3q12o)*0;
    
    dq1spdp1 = dq1spdi*ch*(1+ci) + dq1spdgbysi*(1+ci)*sh - dq1spdhbysi*(1+ci)*sh;
    
    dq1spdp2 = dq1spdi*sh*(1+ci) - dq1spdgbysi*(1+ci)*ch + dq1spdhbysi*(1+ci)*ch;

    dq2spdi = cos(g+h)*(t1isi + t2isi) - sin(g+h)*t3isi;
    
    dq2spdgbysi = (-sin(g+h)*(t1q12o + t2q12o) - cos(g+h)*t3q12o)*0 + (cos(g+h)*(t1g + t2g) - sin(g+h)*t3g);

    dq2spdhbysi = (-sin(g+h)*(t1q12o + t2q12o) - cos(g+h)*t3q12o)*0;
    
    dq2spdp1 = dq2spdi*ch*(1+ci) + dq2spdgbysi*(1+ci)*sh - dq2spdhbysi*(1+ci)*sh;
    
    dq2spdp2 = dq2spdi*sh*(1+ci) - dq2spdgbysi*(1+ci)*ch + dq2spdhbysi*(1+ci)*ch;
    
    % dq1 dq2
    dqdqss1 = a/mu*((e/2*(4*eta^2+9*eta+9)/(1+eta)/eta^4)*cos(f) + 1/2*((7-eta^3-eta^2)/eta^4*cos(2*f)) ...
                    + 3/2*e*cos(3*f)/eta^4 + 1/4*(1-eta^4)/eta^4*cos(4*f) ...
                    + (9 + 9*eta-eta^2-7*eta^3-6*eta^4-4*eta^5)/4/eta^4/(1+eta))*K2n;
    
    dbgterm = getW2k0F(@dqdqss1t2,n,Jn,J2,a,e,eta,i,f,l,G,g,mu,Re);
    dbgt2 = 1/(2*eta^4*L)*(e*(4*eta^2+9*eta+9)*cos(f)/(1+eta) - (eta^3+eta^2-7)*cos(2*f) ...
                + 3*e*cos(3*f) - 1/2*(eta^2 - 1)*cos(4*f) - (6*eta^3+eta^2-9)/2);
       
    [~,~,W2W3fk2bye2q1q1sst2] = getW2W3f(@q1q1sst2nkterm,n,Jn,J2,a,e,eta,i,f,l,G,g,mu,Re);

    q1q1sst2w2k = getW2kF(@q1q1sst2w2kterm,n,1,Jn,J2,a,e,eta,i,f,l,G,g,mu,Re);
            
    dqdqss2 = (2*n-1)/eta^2*W2W3f*dfdl*eta/sqrt(mu*a) + eta/sqrt(mu*a)*d2fdfde*W2W3fk1bye*dfdl ...
                + dbgterm*dbgt2 + W2W3fk2bye2q1q1sst2*eta/L*dfdl ...
                + q1q1sst2w2k*1/2*eta/sqrt(mu*a)*dfdl*(n-1)*sqrt((1-eta)/(1+eta));

    [~,~,W2W3gk2bye2q1q1sst3] = getW2W3g(@q1q1sst3func,n,Jn,J2,a,e,eta,i,f,l,G,g,mu,Re,false);            
    dqdqss3 = -(2*n-1)/G*W2W3g - eta/L*W2W3gk2bye2q1q1sst3;
    
    q1q1sst4 = getW2k0F(@q1q1sst4term,n,Jn,J2,a,e,eta,i,f,l,G,g,mu,Re);
    dqdqss4 = eta/L*dfde*(W2W3f2k1bye*dfdl - W2W3fgk1bye) + eta/L*dfde*dfdlm1bye*q1q1sst4;
    
    [~,q1q1sst5k1bye,~] = getK2ngF(@q1q1sst5term1, n,Jn,J2,a,e,eta,i,g,G,f,mu,Re,false);
    [~,~,q1q1sst5k2bye] = getK2ngF(@q1q1sst5term2, n,Jn,J2,a,e,eta,i,g,G,f,mu,Re,false);
    
    dqdqss5 = Delta*a/mu*(q1q1sst5k1bye + q1q1sst5k2bye);
    
    dqdqssrest = -sqrt(a^3/mu)*dfdlm1bye*K2neGH + 1/eta^2/(1+ci)*(W2iW3f*dfdl - W2iW3g)*eta/L ...
                    - eta/L*sqrt(a^3/mu)*dfde*K2ngbye;
    
    dqdqss = dqdqss1 + dqdqss2 + dqdqss3 + dqdqss4 + dqdqss5 + dqdqssrest;
    
    dqdqcs12t1 = getW2k0F(@q1q1cs12t1,n,Jn,J2,a,e,eta,i,f,l,G,g,mu,Re);
    dqdqcs12t2 = getW2k0F(@q1q1cs12t1,n,Jn,J2,a,e,eta,i,f,l,G,g,mu,Re);
    dqdqcs12t3 = getdqdqcst123(@dummy, n,Jn,J2,a,e,eta,i,g,G,mu,Re);

    dqdqcs1 = (f-l)*a/(mu*eta)*(dqdqcs12t1 + (n-1)/2*dqdqcs12t2 + dqdqcs12t3); 

    dqdqcs2 = a/2/mu/eta^4*(e/(1+eta)*(4*eta^2+7*eta+7)*sin(f) - (eta^3+eta^2-7)*sin(2*f) ...
                + 3*e*sin(3*f) + 0.5*e^2*sin(4*f))*K2n ...
                + K2ngbye*a/(2*mu*eta)*(4*cos(f)+e*cos(2*f) + e*(2*eta^2+3*eta+3)/(1+eta)) ...
                - a*eta/mu*dfdlm1bye*K2ngk1bye;

    dqdqscK3 =   -sqrt(a^3/mu)*(f-l)*K2neGHe + sqrt(a^3/mu)*(-dfde*K2neGH - e/G*dfde*K2n + eta/L*d2fde2*K2n ...
                + dfde*eta/L*K2ne);
    
    dqdqsc4 = -e/G*W2eW3 - e^2/eta^3/L/(1+ci)*W2iW3 -e/G*dfde*W2W3f ...
                + eta/L*(W2e2W3 + 1/eta^2/(1+ci)*W2iW3 + e/eta^2/(1+ci)*W2eiW3 + 2*e^2/eta^4*W2iW3) ...
                + eta/L*(W2eW3f + e/eta^2/(1+ci)*W2iW3f)*dfde ...
                + eta/L*(d2fde2*W2W3f + dfde*W2eW3f + dfde*W2W3f2*dfde);            

    dqdqcs4 = -eta/L/eta^2/(1+ci)*W2iW3;
           
    dfdl2etam1 = 1/eta^5*(e*(eta^4+eta^3+eta^2+eta+1)/(1+eta) + 4*cos(f) + 6*e*(cos(f))^2 ...
                            + 4*e^2*(cos(f))^3 + e^3*(cos(f))^4);
    dqdqcs5sf = 2/L*dfdl2etam1*eta*sin(f)*getW2k0F(@q1q1sc5sfterm,n,Jn,J2,a,e,eta,i,f,l,G,g,mu,Re);
    
    dqdqcs5s1pc = 1/L*getW2k0F(@q1q1sc5s1pcterm,n,Jn,J2,a,e,eta,i,f,l,G,g,mu,Re);
    
    dqdqcs5s1mc = 1/L*getW2k0F(@q1q1sc5s1mcterm,n,Jn,J2,a,e,eta,i,f,l,G,g,mu,Re);
    
    cs5cfgterm = 1/2/eta^4*(e/(1+eta)*(4*eta^2 + 7*eta +7)*sin(f) + (7-eta^2-eta^3)*sin(2*f) ...
                    + 3*e*sin(3*f) + e^2*sin(4*f)/2);
    dqdqcs5cfg = 1/L*cs5cfgterm*getW2k0F(@q1q1sc5cfgterm,n,Jn,J2,a,e,eta,i,f,l,G,g,mu,Re);
    
    dqdqcs5cfg2 = -dfdlm1bye/2/eta/L*(4*cos(f)+e*cos(2*f)+3*e)...
                    *getW2k0F(@q1q1sc5cfg2term,n,Jn,J2,a,e,eta,i,f,l,G,g,mu,Re);
    
    dqdqcs5w3fg = W2W3fgk2bye*eta^2/L*dfdl*(1 + 1/eta);
    dqdqcs5w3g2 = W2W3g2k2bye*eta^2/L*(-1/eta);
    dqdqcs5w3f2 = W2W3f2k2bye*eta^2/L*(-dfdl^2);
    dqdqcs5w3f = -W2W3fk2bye*eta^2/L*d2fdl2;
    dqdqcs5cf = -eta^2/L*(n-1)/2*2*cos(f)*d2fdfdldfdlbye*getW2k0F(@q1q1sc5cfterm,n,Jn,J2,a,e,eta,i,f,l,G,g,mu,Re);
    dqdqcs51pc1mc = -eta^2/L*(n-1)/2*getW2k0F(@q1q1sc51pc1mcterm,n,Jn,J2,a,e,eta,i,f,l,G,g,mu,Re);
    
    dqdqcs5 = -eta/L*(2*n-1)/eta^2*W2W3 + W2W3k2byek + dqdqcs5sf + dqdqcs5s1pc + dqdqcs5s1mc ...
                + (-eta/L)*(dfde*W2W3fk1bye) + dqdqcs5cfg + dqdqcs5cfg2 + dqdqcs5w3fg + dqdqcs5w3g2 + dqdqcs5w3f2 ...
                + dqdqcs5w3f + dqdqcs5cf + dqdqcs51pc1mc;
    
    dqdqcs = dqdqcs1 + dqdqcs2 + dqdqcs4 + dqdqcs5;
    
    dqdqsc = dqdqscK3 + dqdqsc4;
  
    dqdqcc = -W1lgee - W23lgee;
    
    dq1spdq1 = sin(g+h)*sin(g+h)*dqdqss + sin(g+h)*cos(g+h)*dqdqsc ...
                + cos(g+h)*sin(g+h)*dqdqcs +  cos(g+h)*cos(g+h)*dqdqcc;
    
    dq1spdq2 = -sin(g+h)*cos(g+h)*dqdqss + sin(g+h)*sin(g+h)*dqdqsc ...
                - cos(g+h)*cos(g+h)*dqdqcs +  si*cos(g+h)*sin(g+h)*dqdqcc;
            
    dq2spdq1 = -cos(g+h)*sin(g+h)*dqdqss - cos(g+h)*cos(g+h)*dqdqsc ...
                + sin(g+h)*sin(g+h)*dqdqcs +  sin(g+h)*cos(g+h)*dqdqcc;
            
    dq2spdq2 = cos(g+h)*cos(g+h)*dqdqss - cos(g+h)*sin(g+h)*dqdqsc ...
                - sin(g+h)*cos(g+h)*dqdqcs +  sin(g+h)*sin(g+h)*dqdqcc;            
    
            
    % D matrix
    Dspn = [daspda, daspdL, daspdp1, daspdp2, daspdq1, daspdq2;
           dLspda, dLspdL, dLspdp1, dLspdp2, dLspdq1, dLspdq2;
           dp1spda, dp1spdL, dp1spdp1, dp1spdp2, dp1spdq1, dp1spdq2;
           dp2spda, dp2spdL, dp2spdp1, dp2spdp2, dp2spdq1, dp2spdq2;
           dq1spda, dq1spdL, dq1spdp1, dq1spdp2, dq1spdq1, dq1spdq2;
           dq2spda, dq2spdL, dq2spdp1, dq2spdp2, dq2spdq1, dq2spdq2];
       
    Dsp = Dsp + Dspn;
            
end

% % Total First Order SP contributions
% EqSPX1 = J2_SPX1;
% 
% % Total Second Order SP contributions
% EqSPX2 = J2_SPX2 + J2^2/2*Jn_SPX2;
% 
% EqSPX = EqSPX1 + EqSPX2;

EqSP = O2Msign*J2^2/2*EqSP;
Dsp = O2Msign*J2^2/2*Dsp;

end


% function ansp2 = getanSP2(n,Jn,J2,a,e,eta,i,f,l,g,L,G,mu,Re)
% 
% ansp2 = 0;
% if sin(i) == 0
%     return;
% end
% 
% dfdl = (1 + e*cos(f))^2/eta^3;
% K2nSP = getK2nSPF(@dummy, n,Jn,J2,a,e,eta,i,g,G,mu,Re);
% W2W3f = getW2W3f(@dummy,n,Jn,J2,a,e,eta,i,f,l,G,g,mu,Re);
% 
% ansp2 = -2*sqrt(a/mu)*sin(i)*(sqrt(a^3/mu)*(dfdl-1)*K2nSP + W2W3f*dfdl);
% 
% end
% 
% function Lambdasp2 = getLambdaSP2(n,Jn,J2,a,e,eta,i,f,l,g,L,G,mu,Re)
% 
% Lambdasp2 = 0;
% if sin(i) == 0
%     return;
% end
% 
% fLG = -sqrt((1-eta)/(1+eta))*sin(f)*(2+e*cos(f))/G;
% K2nSP = getK2nSPF(@dummy, n,Jn,J2,a,e,eta,i,g,G,mu,Re);
% K2nSPLGH = getK2nSPF(@KLGHfunc, n,Jn,J2,a,e,eta,i,g,G,mu,Re);
% W2W3f = getW2W3f(@dummy,n,Jn,J2,a,e,eta,i,f,l,G,g,mu,Re);
% W2Lt3W3 = getW2W3F(@Lt3func,n,Jn,J2,a,e,eta,i,f,l,G,g,mu,Re);
% W1 = sqrt(a^3/mu)*(f-l)*K2nSP;
% 
% Lambdasp2 = sin(i)*(3/L*W1 + sqrt(a^3/mu)*((f-l)*K2nSPLGH + K2nSP*fLG) - W2Lt3W3 + W2W3f*fLG);
% 
% end
% 
function Lt3 = Lt3func(n,j,k,a,e,eta,i,ci,G,f,g)

Lt3 = 1/G*((2*n-1) + (n-2*j)*ci/(1+ci) + k*eta^2/(1+eta));

end

function Lt3 = Lt3n2jfunc(n,j,k,a,e,eta,i,ci,G,f,g)

Lt3 = Lt3func(n,j,k,a,e,eta,i,G,ci)*(n-2*j);

end

function Lt3 = Lt3i2func(n,j,k,a,e,eta,i,ci,G,f,g)

Lt3 = 1/G*(-sin(i)/(1+ci)^2);

end

function Lt3 = Lt3e2func(n,j,k,a,e,eta,i,ci,G,f,g)

Lt3 = 1/G*(-e*(eta+2)*k/(1+eta^2));

end
% 
% function [p1sp2,p2sp2] = getp1p2SP2(n,Jn,J2,a,e,eta,i,f,l,g,h,L,G,mu,Re)
% 
% p1sp2 = 0; p2sp2 = 0;
% if cos(i) == 0
%     return;
% end
% 
% % K2nSP = getK2nSPF(@dummy, n,Jn,J2,a,eta,i,g,mu,Re);
% K2nSPn2j = getK2nSPF(@n2jfunc, n,Jn,J2,a,e,eta,i,g,G,mu,Re);
% W2n2jW3 = getW2W3F(@n2jfunc,n,Jn,J2,a,e,eta,i,f,l,G,g,mu,Re);
% W2W3g = getW2W3g(@dummy,n,Jn,J2,a,e,eta,i,f,l,G,g,mu,Re);
% K2nSPg = getK2nSPg(n,Jn,J2,a,e,eta,i,g,mu,Re);
% 
% t1 = cos(i)/2/G/((cos(i/2))^2);
% t2 = sqrt(a^3/mu)*(f-l)*(-K2nSPn2j) - W2n2jW3;
% t3 = sqrt(a^3/mu)*(f-l)*K2nSPg + W2W3g;
% 
% p1sp2 = t1*(-sin(h)*t2 - cos(h)*t3);
% p2sp2 = t1*(cos(h)*t2 - sin(h)*t3);
% 
% end
% 
% function [q1sp2,q2sp2] = getq1q2SP2(n,Jn,J2,a,e,eta,i,f,l,g,h,L,G,mu,Re)
% 
% q1sp2 = 0; q2sp2 = 0;
% if sin(i) == 0
%     return;
% end
% 
% dfde = (2+e*cos(f))*sin(f)/eta^2;
% K2nSP = getK2nSPF(@dummy, n,Jn,J2,a,e,eta,i,g,G,mu,Re);
% K2nSPeGH1 = getK2nSPF(@KeGHfunc1, n,Jn,J2,a,e,eta,i,g,G,mu,Re);
% K2nSPeGH2 = getK2nSPFk1(@KeGHfunc2, n,Jn,J2,a,e,eta,i,g,G,mu,Re);
% K2nSPeGH = K2nSPeGH1 + K2nSPeGH2;
% W2eW3 = getW2eW3(n,Jn,J2,a,e,eta,i,f,l,G,g,mu,Re);
% W2iW3SI = getW2iW3SI(n,Jn,J2,a,e,eta,i,f,l,G,g,mu,Re);
% W2W3f = getW2W3f(@dummy,n,Jn,J2,a,e,eta,i,f,l,G,g,mu,Re);
% W1lge = getW1lge(n,Jn,J2,a,e,eta,i,f,l,g,G,mu,Re);
% W23lge = getW23lge(n,Jn,J2,a,e,eta,i,f,l,g,L,G,mu,Re);
% 
% t1 = sqrt(a^3/mu)*((f-l)*(K2nSPeGH) - dfde*eta/L*K2nSP);
% t2 = -eta/L*W2eW3 - e/2/G/((cos(i/2))^2)*W2iW3SI - eta/L*dfde*W2W3f;
% t3 = W1lge + W23lge;
% 
% q1sp2 = sin(i)*(-sin(g+h)*(t1 + t2) - cos(g+h)*t3);
% q2sp2 = sin(i)*(cos(g+h)*(t1 + t2) - sin(g+h)*t3);
% 
% end

function n2j = n2jfunc(n,j,k,a,e,eta,i,ci,G,f,g)

n2j = n-2*j;

end

function k1 = kfunc(n,j,k,a,e,eta,i,ci,G,f,g)

k1 = k;

end

function k1 = kn2jfunc(n,j,k,a,e,eta,i,ci,G,f,g)

k1 = k*(n-2*j);

end

function n2j1 = n2j1func(n,j,k,a,e,eta,i,ci,G,f,g)

n2j1 = n-2*j-1;

end

function term = q1q1sst2nkterm(n,j,k,a,e,eta,i,ci,G,f,g)

term = ((n - k) - eta/k*(n - k))*k/(n-k);

end

function term = q1q1sst3func(n,j,k,a,e,eta,i,ci,G,f,g)

term = (k-1);

end

function n2j1 = n2jn2jfunc(n,j,k,a,e,eta,i,ci,G,f,g)

n2j1 = (n-2*j)^2;

end


function [W2W3,W2W3k1bye,W2W3k2bye] = getW2W3F(func,n,Jn,J2,a,e,eta,i,f,l,G,g,mu,Re,divsi)
ci = cos(i);
% outer j loop
sumj = 0; sumjk1 = 0; sumjk2 = 0;

for j = 0:floor(n/2)
    
    if divsi == true
        beta = getBetaSP(j,n,i);
    else
        beta = getBeta(j,n,i);
    end

    % inner k loop
    sumk = 0; sumkk1 = 0; sumkk2 = 0;
    
    for k = 0:(n-1)
        
        alp = getAlpha(k,n,e);
        
        W3 = getW3(n,j,k,a,e,i,eta,f,g);
        
        sumk = sumk + alp/2^k*W3*func(n,j,k,a,e,eta,i,ci,G,f,g);
        
        if k >= 1
            alp1 = getAlpha1(k,n,e);
            sumkk1 = sumkk1 + (n-k)/k*alp1/2^k*W3*func(n,j,k,a,e,eta,i,ci,G,f,g);
        end
        
        if k >= 2
            alp2 = getAlpha2(k,n,e);
            sumkk2 = sumkk2 + (n-k)/k*alp2/2^k*W3*func(n,j,k,a,e,eta,i,ci,G,f,g);
        end        
    end
    
    sumj = sumj + beta*sumk;
    sumjk1 = sumjk1 + beta*sumkk1;
    sumjk2 = sumjk2 + beta*sumkk2;

end

del = getDeltad(n,Jn,J2,a,eta,mu,Re);
W2W3 = del*sumj;
W2W3k1bye = del*sumjk1;
W2W3k2bye = del*sumjk2;
end

function W2W3F = getW2W3m1F(func,n,Jn,J2,a,e,eta,i,f,l,G,g,mu,Re)

% outer j loop
sumj = 0;
for j = 0:floor((n-2)/2)
    betaSP = getBetaSPm1(j,n,i);

    % inner k loop
    sumk = 0;
    for k = 0:(n-1)
        
        alp = getAlpha(k,n,e);
        W3 = getW3(n,j,k,a,e,i,eta,f,g);
        sumk = sumk + alp/2^k*W3*func(n,j,k,a,e,eta,i,ci,G,f,g);
    end
    if isinf(betaSP) == 0
        sumj = sumj + betaSP*sumk;
    elseif sumk ~= 0
        disp('Inclination singularity');
    end        
    
end

del = getDeltad(n,Jn,J2,a,eta,mu,Re);
W2W3F = del*sumj;
end


function k0dadldadg = k0dadldadgterm(n,j,k,a,e,eta,i,ci,G,f,g)

odd = 0;
if mod(n,2) ~= 0
    % odd harmonic
    odd = 1;
end

sump = 0;

for p = 0:(floor((n-1)/2)-j)
    
    c = n-2*j-2*p;
    
    gam = getGamma(p,j,n);
    
    if odd ~= 1
        sump = sump - gam*sin(c*f+c*g)*c;
    else
        sump = sump + gam*cos(c*f+c*g)*c;
    end
end

k0dadldadg = sump;

end


function q1q1sst4t = q1q1sst4term(n,j,k,a,e,eta,i,ci,G,f,g)

odd = 0;
if mod(n,2) ~= 0
    % odd harmonic
    odd = 1;
end

sump = 0;

for p = 0:(floor((n-1)/2)-j)
    
    c = n-2*j-2*p;
    
    gam = getGamma(p,j,n);
    
    if odd ~= 1
        sump = sump - gam*sin(c*f+c*g)*c;
    else
        sump = sump + gam*cos(c*f+c*g)*c;
    end
end

q1q1sst4t = sump;

end

function k0W2Lt3W3fg = k0W2Lt3W3fgterm(n,j,k,a,e,eta,i,ci,G,f,g)

odd = 0;
if mod(n,2) ~= 0
    % odd harmonic
    odd = 1;
end

sump = 0;

Lt3 = Lt3func(n,j,k,a,e,eta,i,G,ci);

for p = 0:(floor((n-1)/2)-j)
    c = n-2*j-2*p;
    gam = getGamma(p,j,n);
    if odd ~= 1
        sump = sump + gam*cos(c*f+c*g)*c;
    else
        sump = sump + gam*sin(c*f+c*g)*c;
    end
end

k0W2Lt3W3fg = sump*Lt3;

end

function dp1spdlgk0term = dp1spdlgk0func(n,j,k,a,e,eta,i,ci,G,f,g)

odd = 0;
if mod(n,2) ~= 0
    % odd harmonic
    odd = 1;
end

sump = 0;

for p = 0:(floor((n-1)/2)-j)
    c = n-2*j-2*p;
    gam = getGamma(p,j,n);
    if odd ~= 1
        sump = sump + gam*cos(c*f+c*g);
    else
        sump = sump + gam*sin(c*f+c*g);
    end
end

dp1spdlgk0term = sump*(n-2*j);

end

function dp1spdlg2k0term = dp1spdlg2k0func(n,j,k,a,e,eta,i,ci,G,f,g)

odd = 0;
if mod(n,2) ~= 0
    % odd harmonic
    odd = 1;
end

sump = 0;

for p = 0:(floor((n-1)/2)-j)
    c = n-2*j-2*p;
    gam = getGamma(p,j,n);
    if odd ~= 1
        sump = sump + gam*sin(c*f+c*g)*c;
    else
        sump = sump - gam*cos(c*f+c*g)*c;
    end
end

dp1spdlg2k0term = sump;

end

function W2k0F = getW2k0F(func,n,Jn,J2,a,e,eta,i,f,l,G,g,mu,Re)

ci = cos(i);

% outer j loop
sumj = 0;
for j = 0:floor(n/2)
    
    beta = getBeta(j,n,i);

    sumk = func(n,j,0,a,e,eta,i,ci,G,f,g);

    sumj = sumj + beta*sumk;
    
end

del = getDeltad(n,Jn,J2,a,eta,mu,Re);
W2k0F = del*sumj;

end

function W2kF = getW2kF(func,n,k,Jn,J2,a,e,eta,i,f,l,G,g,mu,Re)

ci = cos(i);

% outer j loop
sumj = 0;
for j = 0:floor(n/2)
    
    beta = getBeta(j,n,i);

    sumk = func(n,j,k,a,e,eta,i,ci,G,f,g);

    sumj = sumj + beta*sumk;
    
end

del = getDeltad(n,Jn,J2,a,eta,mu,Re);
W2kF = del*sumj;

end


function [W2W3fg,W2W3fgk1bye,W2W3fgk2bye] = getW2W3fg(func,n,Jn,J2,a,e,eta,i,f,l,G,g,mu,Re,divsi)

ci = cos(i);

% outer j loop
sumj = 0; sumjk1 = 0; sumjk2 = 0;

for j = 0:floor(n/2)
    
    if divsi == true
        beta = getBetaSP(j,n,i);
    else
        beta = getBeta(j,n,i);
    end

    % inner k loop
    sumk = 0; sumkk1 = 0; sumkk2 = 0;
    
    for k = 0:(n-1)
        
        alp = getAlpha(k,n,e);
        
        W3fg = getW3fg(n,j,k,a,e,i,eta,f,g);
        
        sumk = sumk + alp/2^k*W3fg*func(n,j,k,a,e,eta,i,ci,G,f,g);
        
        if k >= 1
            alp1 = getAlpha1(k,n,e);
            sumkk1 = sumkk1 + (n-k)/k*alp1/2^k*W3fg*func(n,j,k,a,e,eta,i,ci,G,f,g);
        end
        
        if k >= 2
            alp2 = getAlpha2(k,n,e);
            sumkk2 = sumkk2 + (n-k)/k*alp2/2^k*W3fg*func(n,j,k,a,e,eta,i,ci,G,f,g);
        end        
    end
    
    sumj = sumj + beta*sumk;
    sumjk1 = sumjk1 + beta*sumkk1;
    sumjk2 = sumjk2 + beta*sumkk2;

end

del = getDeltad(n,Jn,J2,a,eta,mu,Re);
W2W3fg = del*sumj;
W2W3fgk1bye = del*sumjk1;
W2W3fgk2bye = del*sumjk2;
end

function [W2W3f,W2W3fk1bye,W2W3fk2bye] = getW2W3f(func,n,Jn,J2,a,e,eta,i,f,l,G,g,mu,Re)
ci = cos(i);
% outer j loop
sumj = 0; sumjk1 = 0; sumjk2 = 0;

for j = 0:floor(n/2)
    
    beta = getBeta(j,n,i);

    % inner k loop
    sumk = 0; sumkk1 = 0; sumkk2 = 0;
    
    for k = 0:(n-1)
        
        alp = getAlpha(k,n,e);
        
        W3f = getW3f(n,j,k,a,e,i,eta,f,g);
        
        sumk = sumk + alp/2^k*W3f*func(n,j,k,a,e,eta,i,ci,G,f,g);
        
        if k >= 1
            alp1 = getAlpha1(k,n,e);
            sumkk1 = sumkk1 + (n-k)/k*alp1/2^k*W3f*func(n,j,k,a,e,eta,i,ci,G,f,g);
        end
        
        if k >= 2
            alp2 = getAlpha2(k,n,e);
            sumkk2 = sumkk2 + (n-k)/k*alp2/2^k*W3f*func(n,j,k,a,e,eta,i,ci,G,f,g);
        end        
    end
    
    sumj = sumj + beta*sumk;
    sumjk1 = sumjk1 + beta*sumkk1;
    sumjk2 = sumjk2 + beta*sumkk2;

end

del = getDeltad(n,Jn,J2,a,eta,mu,Re);
W2W3f = del*sumj;
W2W3fk1bye = del*sumjk1;
W2W3fk2bye = del*sumjk2;
end


function [W2W3f2,W2W3f2k1bye,W2W3f2k2bye] = getW2W3f2(func,n,Jn,J2,a,e,eta,i,f,l,G,g,mu,Re,divsi)

ci = cos(i);

% outer j loop
sumj = 0; sumjk1 = 0; sumjk2 = 0;

for j = 0:floor(n/2)
    
    if divsi == true
        beta = getBetaSP(j,n,i);
    else
        beta = getBeta(j,n,i);
    end

    % inner k loop
    sumk = 0; sumkk1 = 0; sumkk2 = 0;
    
    for k = 0:(n-1)
        
        alp = getAlpha(k,n,e);
        
        W3f2 = getW3f2(n,j,k,a,e,i,eta,f,g);
        
        sumk = sumk + alp/2^k*W3f2*func(n,j,k,a,e,eta,i,ci,G,f,g);
        
        if k >= 1
            alp1 = getAlpha1(k,n,e);
            sumkk1 = sumkk1 + (n-k)/k*alp1/2^k*W3f2*func(n,j,k,a,e,eta,i,ci,G,f,g);
        end
        
        if k >= 2
            alp2 = getAlpha2(k,n,e);
            sumkk2 = sumkk2 + (n-k)/k*alp2/2^k*W3f2*func(n,j,k,a,e,eta,i,ci,G,f,g);
        end        
    end
    
    sumj = sumj + beta*sumk;
    sumjk1 = sumjk1 + beta*sumkk1;
    sumjk2 = sumjk2 + beta*sumkk2;

end

del = getDeltad(n,Jn,J2,a,eta,mu,Re);
W2W3f2 = del*sumj;
W2W3f2k1bye = del*sumjk1;
W2W3f2k2bye = del*sumjk2;
end


function [W2W3g,W2W3gk1bye,W2W3gk2bye] = getW2W3g(func,n,Jn,J2,a,e,eta,i,f,l,G,g,mu,Re,divsi)
ci = cos(i);
% outer j loop
sumj = 0; sumjk1 = 0; sumjk2 = 0;

for j = 0:floor(n/2)
    
    if divsi == true
        beta = getBetaSP(j,n,i);
    else
        beta = getBeta(j,n,i);
    end

    % inner k loop
    sumk = 0; sumkk1 = 0; sumkk2 = 0;
    
    for k = 0:(n-1)
        
        alp = getAlpha(k,n,e);
        
        W3g = getW3g(n,j,k,a,e,i,eta,f,g);
        
        sumk = sumk + alp/2^k*W3g*func(n,j,k,a,e,eta,i,ci,G,f,g);
        
        if k >= 1
            alp1 = getAlpha1(k,n,e);
            sumkk1 = sumkk1 + (n-k)/k*alp1/2^k*W3g*func(n,j,k,a,e,eta,i,ci,G,f,g);
        end
        
        if k >= 2
            alp2 = getAlpha2(k,n,e);
            sumkk2 = sumkk2 + (n-k)/k*alp2/2^k*W3g*func(n,j,k,a,e,eta,i,ci,G,f,g);
        end        
    end
    
    sumj = sumj + beta*sumk;
    sumjk1 = sumjk1 + beta*sumkk1;
    sumjk2 = sumjk2 + beta*sumkk2;

end

del = getDeltad(n,Jn,J2,a,eta,mu,Re);
W2W3g = del*sumj;
W2W3gk1bye = del*sumjk1;
W2W3gk2bye = del*sumjk2;
end

function [W2W3g2,W2W3g2k1bye,W2W3g2k2bye] = getW2W3g2(func,n,Jn,J2,a,e,eta,i,f,l,G,g,mu,Re,divsi)
ci = cos(i);
% outer j loop
sumj = 0; sumjk1 = 0; sumjk2 = 0;

for j = 0:floor(n/2)
    
    if divsi == true
        beta = getBetaSP(j,n,i);
    else
        beta = getBeta(j,n,i);
    end

    % inner k loop
    sumk = 0; sumkk1 = 0; sumkk2 = 0;
    
    for k = 0:(n-1)
        
        alp = getAlpha(k,n,e);
        
        W3g2 = getW3g2(n,j,k,a,e,i,eta,f,g);
        
        sumk = sumk + alp/2^k*W3g2*func(n,j,k,a,e,eta,i,ci,G,f,g);
        
        if k >= 1
            alp1 = getAlpha1(k,n,e);
            sumkk1 = sumkk1 + (n-k)/k*alp1/2^k*W3g2*func(n,j,k,a,e,eta,i,ci,G,f,g);
        end
        
        if k >= 2
            alp2 = getAlpha2(k,n,e);
            sumkk2 = sumkk2 + (n-k)/k*alp2/2^k*W3g2*func(n,j,k,a,e,eta,i,ci,G,f,g);
        end        
    end
    
    sumj = sumj + beta*sumk;
    sumjk1 = sumjk1 + beta*sumkk1;
    sumjk2 = sumjk2 + beta*sumkk2;

end

del = getDeltad(n,Jn,J2,a,eta,mu,Re);

W2W3g2 = del*sumj;
W2W3g2k1bye = del*sumjk1;
W2W3g2k2bye = del*sumjk2;

end



function [W2W3g,W2W3gk1bye] = getW2W3gm1F(func,n,Jn,J2,a,e,eta,i,f,l,G,g,mu,Re)

% outer j loop
sumj = 0; sumjk1 = 0;
for j = 0:floor((n-2)/2)
    betaSP = getBetaSPm1(j,n,i);

    % inner k loop
    sumk = 0; sumkk1 = 0;
    for k = 0:(n-1)
        
        alp = getAlpha(k,n,e);
        W3g = getW3g(n,j,k,a,e,i,eta,f,g);
        sumk = sumk + alp/2^k*W3g*func(n,j,k,a,e,eta,i,ci,G,f,g);
        
        if k >= 1
            alp1 = getAlpha1(k,n,e);
            sumkk1 = sumkk1 + (n-k)/k*alp1/2^k*W3g*func(n,j,k,a,e,eta,i,ci,G,f,g);
        end
        
    end
    if isinf(betaSP) == 0
        sumj = sumj + betaSP*sumk;
        sumjk1 = sumjk1 + betaSP*sumkk1;
    elseif sumk ~= 0
        disp('Inclination singularity');
    end        
    
end

del = getDeltad(n,Jn,J2,a,eta,mu,Re);
W2W3g = del*sumj;
W2W3gk1bye = del*sumjk1;
end

function [W2W3g,W2W3gk1bye] = getW2W3g2m1F(func,n,Jn,J2,a,e,eta,i,f,l,G,g,mu,Re)

% outer j loop
sumj = 0; sumjk1 = 0;
for j = 0:floor((n-2)/2)
    betaSP = getBetaSPm1(j,n,i);

    % inner k loop
    sumk = 0; sumkk1 = 0;
    for k = 0:(n-1)
        
        alp = getAlpha(k,n,e);
        W3g2 = getW3g2(n,j,k,a,e,i,eta,f,g);
        sumk = sumk + alp/2^k*W3g2*func(n,j,k,a,e,eta,i,ci,G,f,g);
        
        if k >= 1
            alp1 = getAlpha1(k,n,e);
            sumkk1 = sumkk1 + (n-k)/k*alp1/2^k*W3g2*func(n,j,k,a,e,eta,i,ci,G,f,g);
        end
        
    end
    if isinf(betaSP) == 0
        sumj = sumj + betaSP*sumk;
        sumjk1 = sumjk1 + betaSP*sumkk1;
    elseif sumk ~= 0
        disp('Inclination singularity');
    end        
    
end

del = getDeltad(n,Jn,J2,a,eta,mu,Re);
W2W3g = del*sumj;
W2W3gk1bye = del*sumjk1;
end

function W2eW3 = getW2eW3F(func,n,Jn,J2,a,e,eta,i,f,l,G,g,mu,Re,divsi)

ci = cos(i);
% outer j loop
sumj = 0;
for j = 0:floor(n/2)
    
    if divsi == true
        beta = getBetaSP(j,n,i);
    else
        beta = getBeta(j,n,i);
    end

    % inner k loop
    sumk1 = 0;
    sumk2 = 0;
    for k = 0:(n-1)
        
        alp = getAlpha(k,n,e);
        W3 = getW3(n,j,k,a,e,i,eta,f,g);
        
        sumk1 = sumk1 + alp/2^k*W3*func(n,j,k,a,e,eta,i,ci,G,f,g);
        
        if k > 0
            alp1 = getAlpha1(k,n,e);
            sumk2 = sumk2 + alp1/2^k*(n-k)*W3*func(n,j,k,a,e,eta,i,ci,G,f,g);
        end
    end
    sumj = sumj + beta*((2*n-1)*e/eta^2*sumk1 + sumk2);

end

del = getDeltad(n,Jn,J2,a,eta,mu,Re);
W2eW3 = del*sumj;

end

function W2eW3f = getW2eW3fF(func,n,Jn,J2,a,e,eta,i,f,l,G,g,mu,Re,divsi)

ci = cos(i);

% outer j loop
sumj = 0;
for j = 0:floor(n/2)
    
    if divsi == true
        beta = getBetaSP(j,n,i);
    else
        beta = getBeta(j,n,i);
    end

    % inner k loop
    sumk1 = 0;
    sumk2 = 0;
    for k = 0:(n-1)
        
        alp = getAlpha(k,n,e);
        W3f = getW3f(n,j,k,a,e,i,eta,f,g);
        
        sumk1 = sumk1 + alp/2^k*W3f*func(n,j,k,a,e,eta,i,ci,G,f,g);
        
        if k > 0
            alp1 = getAlpha1(k,n,e);
            sumk2 = sumk2 + alp1/2^k*(n-k)*W3f*func(n,j,k,a,e,eta,i,ci,G,f,g);
        end
    end
    sumj = sumj + beta*((2*n-1)*e/eta^2*sumk1 + sumk2);

end

del = getDeltad(n,Jn,J2,a,eta,mu,Re);
W2eW3f = del*sumj;

end

function W2eW3g = getW2eW3gF(func,n,Jn,J2,a,e,eta,i,f,l,G,g,mu,Re,divsi)

ci = cos(i);

% outer j loop
sumj = 0;
for j = 0:floor(n/2)
    
    if divsi == true
        beta = getBetaSP(j,n,i);
    else
        beta = getBeta(j,n,i);
    end

    % inner k loop
    sumk1 = 0;
    sumk2 = 0;
    for k = 0:(n-1)
        
        alp = getAlpha(k,n,e);
        W3g = getW3g(n,j,k,a,e,i,eta,f,g);
        
        sumk1 = sumk1 + alp/2^k*W3g*func(n,j,k,a,e,eta,i,ci,G,f,g);
        
        if k > 0
            alp1 = getAlpha1(k,n,e);
            sumk2 = sumk2 + alp1/2^k*(n-k)*W3g*func(n,j,k,a,e,eta,i,ci,G,f,g);
        end
    end
    sumj = sumj + beta*((2*n-1)*e/eta^2*sumk1 + sumk2);

end

del = getDeltad(n,Jn,J2,a,eta,mu,Re);
W2eW3g = del*sumj;

end



% function W2eW3F = getW2eW3F(func,n,Jn,J2,a,e,eta,i,f,l,G,g,mu,Re)
% 
% % outer j loop
% sumj = 0;
% for j = 0:floor(n/2)
%     betaSP = getBetaSP(j,n,i);
% 
%     % inner k loop
%     sumk1 = 0;
%     sumk2 = 0;
%     for k = 0:(n-1)
%         
%         alp = getAlpha(k,n,e);
%         W3 = getW3(n,j,k,a,e,i,eta,f,g);
%         sumk1 = sumk1 + alp/2^k*W3f;
%         
%         if k > 0
%             alp1 = getAlpha1(k,n,e);
%             sumk2 = sumk2 + alp1/2^k*(n-k)*W3*func(n,j,k,a,e,eta,i,ci,G,f,g);
%         end
%     end
%     if isinf(betaSP) == 0
%         sumj = sumj + betaSP*((2*n-1)*e/eta^2*sumk1 + sumk2);
%     elseif sumk1 ~= 0 && sumk2 ~= 0
%         disp('Inclination singularity');
%     end
% 
% end
% 
% del = getDeltad(n,Jn,J2,a,eta,mu,Re);
% W2eW3F = del*sumj;
% end

function W2eW3F = getW2e2W3F(func,n,Jn,J2,a,e,eta,i,f,l,G,g,mu,Re)

ci = cos(i);

% outer j loop
sumj = 0;
for j = 0:floor(n/2)
    beta = getBeta(j,n,i);

    % inner k loop
    sumk1 = 0;
    sumk2 = 0;
    sumk3 = 0;
    for k = 0:(n-1)
        
        alp = getAlpha(k,n,e);
        W3 = getW3(n,j,k,a,e,i,eta,f,g);
        sumk1 = sumk1 + alp/2^k*W3;
        
        if k > 0
            alp1 = getAlpha1(k,n,e);
            sumk2 = sumk2 + alp1/2^k*(n-k)*W3*func(n,j,k,a,e,eta,i,ci,G,f,g);
        end
        if k > 1
            alp2 = getAlpha2(k,n,e);
            sumk3 = sumk3 + alp2/2^k*(n-k)*W3*func(n,j,k,a,e,eta,i,ci,G,f,g);
        end        
    end
    sumj = sumj + beta*((2*n-1)*(1+2*n*e^2)/eta^4*sumk1 + 2*(2*n-1)*e/eta^2*sumk2 + sumk3);
end

del = getDeltad(n,Jn,J2,a,eta,mu,Re);
W2eW3F = del*sumj;
end


function W2iW3 = getW2iW3(func,n,Jn,J2,a,e,eta,i,f,l,G,g,mu,Re,divsi)

ci = cos(i);

% outer j loop
sumj = 0;
for j = 0:floor(n/2)
    
    if divsi == true
        beta = getBetaSP(j,n,i);
    else
        beta = getBeta(j,n,i);
    end

    % inner k loop
    sumk1 = 0;
    for k = 0:(n-1)
        
        alp = getAlpha(k,n,e);
        W3 = getW3(n,j,k,a,e,i,eta,f,g);
        sumk1 = sumk1 + alp/2^k*W3*func(n,j,k,a,e,eta,i,ci,G,f,g);
        
    end

    sumj = sumj + (n-2*j)*beta*sumk1;

end

del = getDeltad(n,Jn,J2,a,eta,mu,Re);
W2iW3 = cos(i)*del*sumj;
end

function W2iW3f = getW2iW3f(func,n,Jn,J2,a,e,eta,i,f,l,G,g,mu,Re,divsi)

ci = cos(i);

% outer j loop
sumj = 0;
for j = 0:floor(n/2)
    
    if divsi == true
        beta = getBetaSP(j,n,i);
    else
        beta = getBeta(j,n,i);
    end

    % inner k loop
    sumk1 = 0;
    for k = 0:(n-1)
        
        alp = getAlpha(k,n,e);
        W3f = getW3f(n,j,k,a,e,i,eta,f,g);
        sumk1 = sumk1 + alp/2^k*W3f*func(n,j,k,a,e,eta,i,ci,G,f,g);
        
    end

    sumj = sumj + (n-2*j)*beta*sumk1;

end

del = getDeltad(n,Jn,J2,a,eta,mu,Re);
W2iW3f = cos(i)*del*sumj;
end


function W2iW3g = getW2iW3g(func,n,Jn,J2,a,e,eta,i,f,l,G,g,mu,Re,divsi)

ci = cos(i);

% outer j loop
sumj = 0;
for j = 0:floor(n/2)
    
    if divsi == true
        beta = getBetaSP(j,n,i);
    else
        beta = getBeta(j,n,i);
    end

    % inner k loop
    sumk1 = 0;
    for k = 0:(n-1)
        
        alp = getAlpha(k,n,e);
        W3g = getW3g(n,j,k,a,e,i,eta,f,g);
        sumk1 = sumk1 + alp/2^k*W3g*func(n,j,k,a,e,eta,i,ci,G,f,g);
        
    end

    sumj = sumj + (n-2*j)*beta*sumk1;

end

del = getDeltad(n,Jn,J2,a,eta,mu,Re);
W2iW3g = cos(i)*del*sumj;
end

function W2i2W3SI2 = getW2i2W3SI2(func,n,Jn,J2,a,e,eta,i,f,l,G,g,mu,Re)

% outer j loop
sumj = 0;
for j = 0:floor(n/2)
    betaSP = getBetaSP(j,n,i);

    % inner k loop
    sumk1 = 0;
    for k = 0:(n-1)
        
        alp = getAlpha(k,n,e);
        W3 = getW3(n,j,k,a,e,i,eta,f,g);
        sumk1 = sumk1 + alp/2^k*W3*func(n,j,k,a,e,eta,i,ci,G,f,g);
        
    end
    if isinf(betaSP) == 0
        sumj = sumj + (n-2*j)*((n-2*j-1)*(cos(i))^2-1)*betaSP*sumk1;
    elseif sumk1 ~= 0
        disp('Inclination singularity');
    end
    

end

del = getDeltad(n,Jn,J2,a,eta,mu,Re);
W2i2W3SI2 = cos(i)*del*sumj;
end


function [dbsi] = delbetasiF(func,n,k,Jn,J2,a,e,eta,i,f,l,G,g,mu,Re)

% outer j loop
sumj = 0;

for j = 0:floor(n/2)
    beta = getBeta(j,n,i);

    sumj = sumj + beta*func(n,j,k,a,e,eta,i,ci,G,f,g);
    
end

del = getDeltad(n,Jn,J2,a,eta,mu,Re);
dbsi = del*sumj;

end


function [dbgsi] = delbetagammasiF(func,n,k,Jn,J2,a,e,eta,i,f,l,G,g,mu,Re)

odd = 0;
if mod(n,2) ~= 0
    % odd harmonic
    odd = 1;
end


% outer j loop
sumj = 0;

for j = 0:floor(n/2)
    beta = getBeta(j,n,i);

    % inner p loop
    sump = 0; 
    for p = 0:(floor(((n-1)/2)-j))
        
        gam = getGamma(p,j,n);
        c = n - 2*j - 2*p;
        if odd == 0 
            sump = sump + gam*cos(c*f + c*g)*func(n,j,k,a,e,eta,i,ci,G,f,g);
        else
            sump = sump + gam*sin(c*f + c*g)*func(n,j,k,a,e,eta,i,ci,G,f,g);
        end
    end

    sumj = sumj + beta*sump;
    
end

del = getDeltad(n,Jn,J2,a,eta,mu,Re);
dbgsi = del*sumj;

end

function [dbg] = dqdqss1t2(n,j,k,a,e,eta,i,ci,G,f,g)

odd = 0;
if mod(n,2) ~= 0
    % odd harmonic
    odd = 1;
end

% inner p loop
sump = 0;
for p = 0:(floor(((n-1)/2)-j))
    
    gam = getGamma(p,j,n);
    c = n - 2*j - 2*p;
    if odd == 0
        sump = sump + gam*cos(c*f + c*g);
    else
        sump = sump + gam*sin(c*f + c*g);
    end
end

dbg = sump;
end


function [dbgsi] = delbetagammaSCsiF(func,n,Jn,J2,a,e,eta,i,f,l,G,g,mu,Re)

odd = 0;
if mod(n,2) ~= 0
    % odd harmonic
    odd = 1;
end


% outer j loop
sumj = 0;

for j = 0:floor(n/2)
    beta = getBeta(j,n,i);

    % inner p loop
    sump = 0; 
    for p = 0:(floor(((n-1)/2)-j))
        
        gam = getGamma(p,j,n);
        c = n - 2*j - 2*p;
        if odd == 0 
            sump = sump - gam*sin(c*f + c*g)*c*func(n,j,k,a,e,eta,i,ci,G,f,g);
        else
            sump = sump + gam*cos(c*f + c*g)*c*func(n,j,k,a,e,eta,i,ci,G,f,g);
        end
    end

    sumj = sumj + beta*sump;
    
end

del = getDeltad(n,Jn,J2,a,eta,mu,Re);
dbgsi = del*sumj;

end

function W3 = getW3(n,j,k,a,e,i,eta,f,g)

odd = 0;
if mod(n,2) ~= 0
    % odd harmonic
    odd = 1;
end

lk = 0;
if mod(k,2) == 0
    lk = 1;
end

t1 = 0;
if odd == 0
    C = getC(k,f);
    t1 = nchoosek(n-2*j, n/2-j)*C;
end

sump = 0;
for p = 0:(floor(((n-1)/2)-j))

gam = getGamma(p,j,n);

sums1 = 0;
sums2 = 0;

c = n - 2*j - 2*p;
sstar = (k - c)/2;
    
for s = 0:floor((k-1)/2)
   
    bcks = nchoosek(k,s);
    
    if odd == 0
        sums1 = sums1 + bcks*(sin((k-2*s+c)*f + c*g))/(k-2*s+c);
        if s ~= sstar
            sums2 = sums2 + bcks*(sin((k-2*s-c)*f - c*g))/(k-2*s-c);
        end
    else
        sums1 = sums1 - bcks*(cos((k-2*s+c)*f + c*g))/(k-2*s+c);
        if s ~= sstar
            sums2 = sums2 + bcks*(cos((k-2*s-c)*f - c*g))/(k-2*s-c);
        end
    end
    
end

t23 = 0;
if lk == 1
    if odd == 0
        t23 = nchoosek(k,k/2)*sin(c*f + c*g)/c;
    else
        t23 = -nchoosek(k,k/2)*cos(c*f + c*g)/c;
    end
end

sump = sump + gam*(sums1 + sums2 + t23);

end

W3 = t1 + sump;

end

function W3f = getW3f(n,j,k,a,e,i,eta,f,g)

odd = 0;
if mod(n,2) ~= 0
    % odd harmonic
    odd = 1;
end

lk = 0;
if mod(k,2) == 0
    lk = 1;
end

t1 = 0;
if odd == 0
    Cf = getCf(k,f);
    t1 = nchoosek(n-2*j, n/2-j)*Cf;
end

sump = 0;
for p = 0:(floor(((n-1)/2)-j))

gam = getGamma(p,j,n);

sums1 = 0;
sums2 = 0;

c = n - 2*j - 2*p;
sstar = (k - c)/2;
    
for s = 0:floor((k-1)/2)
   
    bcks = nchoosek(k,s);

    if odd == 0
        sums1 = sums1 + bcks*(cos((k-2*s+c)*f + c*g));
        if s ~= sstar
            sums2 = sums2 + bcks*(cos((k-2*s-c)*f - c*g));
        end
    else
        sums1 = sums1 + bcks*(sin((k-2*s+c)*f + c*g));
        if s ~= sstar
            sums2 = sums2 - bcks*(sin((k-2*s-c)*f - c*g));
        end
    end
    
end

t23 = 0;
if lk == 1
    if odd == 0
        t23 = nchoosek(k,k/2)*cos(c*f + c*g);
    else
        t23 = nchoosek(k,k/2)*sin(c*f + c*g);
    end
end

sump = sump + gam*(sums1 + sums2 + t23);

end

W3f = t1 + sump;

end

function W3f2 = getW3f2(n,j,k,a,e,i,eta,f,g)

odd = 0;
if mod(n,2) ~= 0
    % odd harmonic
    odd = 1;
end

lk = 0;
if mod(k,2) == 0
    lk = 1;
end

t1 = 0;
if odd == 0
    Sf = getSf(k,f);
    t1 = nchoosek(n-2*j, n/2-j)*Sf;
end

sump = 0;
for p = 0:(floor(((n-1)/2)-j))

gam = getGamma(p,j,n);

sums1 = 0;
sums2 = 0;

c = n - 2*j - 2*p;
sstar = (k - c)/2;
    
for s = 0:floor((k-1)/2)
   
    bcks = nchoosek(k,s);

    if odd == 0
        sums1 = sums1 - bcks*(sin((k-2*s+c)*f + c*g))*(k-2*s+c);
        if s ~= sstar
            sums2 = sums2 - bcks*(sin((k-2*s-c)*f - c*g))*(k-2*s-c);
        end
    else
        sums1 = sums1 + bcks*(cos((k-2*s+c)*f + c*g))*(k-2*s+c);
        if s ~= sstar
            sums2 = sums2 - bcks*(cos((k-2*s-c)*f - c*g))*(k-2*s-c);
        end
    end
    
end

t23 = 0;
if lk == 1
    if odd == 0
        t23 = -nchoosek(k,k/2)*sin(c*f + c*g)*c;
    else
        t23 = nchoosek(k,k/2)*cos(c*f + c*g)*c;
    end
end

sump = sump + gam*(sums1 + sums2 + t23);

end

W3f2 = t1 + sump;

end

function W3g = getW3g(n,j,k,a,e,i,eta,f,g)

odd = 0;
if mod(n,2) ~= 0
    % odd harmonic
    odd = 1;
end

lk = 0;
if mod(k,2) == 0
    lk = 1;
end

sump = 0;
for p = 0:(floor(((n-1)/2)-j))

gam = getGamma(p,j,n);

sums1 = 0;
sums2 = 0;

c = n - 2*j - 2*p;
sstar = (k - c)/2;
    
for s = 0:floor((k-1)/2)

    bcks = nchoosek(k,s);

    if odd == 0
        sums1 = sums1 + bcks*(cos((k-2*s+c)*f + c*g))*c/(k-2*s+c);
        if s ~= sstar
            sums2 = sums2 - bcks*(cos((k-2*s-c)*f - c*g))*c/(k-2*s-c);
        end
    else
        sums1 = sums1 + bcks*(sin((k-2*s+c)*f + c*g))*c/(k-2*s+c);
        if s ~= sstar
            sums2 = sums2 + bcks*(sin((k-2*s-c)*f - c*g))*c/(k-2*s-c);
        end
    end
    
end

t23 = 0;
if lk == 1
    if odd == 0
        t23 = nchoosek(k,k/2)*cos(c*f + c*g);
    else
        t23 = nchoosek(k,k/2)*sin(c*f + c*g);
    end
end

sump = sump + gam*(sums1 + sums2 + t23);

end

W3g = sump;

end

function W3g = getW3g2(n,j,k,a,e,i,eta,f,g)

odd = 0;
if mod(n,2) ~= 0
    % odd harmonic
    odd = 1;
end

lk = 0;
if mod(k,2) == 0
    lk = 1;
end

sump = 0;
for p = 0:(floor(((n-1)/2)-j))

gam = getGamma(p,j,n);

sums1 = 0;
sums2 = 0;

c = n - 2*j - 2*p;
sstar = (k - c)/2;
    
for s = 0:floor((k-1)/2)

    bcks = nchoosek(k,s);

    if odd == 0
        sums1 = sums1 - bcks*(sin((k-2*s+c)*f + c*g))*c^2/(k-2*s+c);
        if s ~= sstar
            sums2 = sums2 + bcks*(sin((k-2*s-c)*f - c*g))*c^2/(k-2*s-c);
        end
    else
        sums1 = sums1 + bcks*(cos((k-2*s+c)*f + c*g))*c^2/(k-2*s+c);
        if s ~= sstar
            sums2 = sums2 + bcks*(cos((k-2*s-c)*f - c*g))*c^2/(k-2*s-c);
        end
    end
    
end

t23 = 0;
if lk == 1
    if odd == 0
        t23 = -nchoosek(k,k/2)*sin(c*f + c*g)*c;
    else
        t23 = nchoosek(k,k/2)*cos(c*f + c*g)*c;
    end
end

sump = sump + gam*(sums1 + sums2 + t23);

end

W3g = sump;

end


function W3fg = getW3fg(n,j,k,a,e,i,eta,f,g)

odd = 0;
if mod(n,2) ~= 0
    % odd harmonic
    odd = 1;
end

lk = 0;
if mod(k,2) == 0
    lk = 1;
end

sump = 0;
for p = 0:(floor(((n-1)/2)-j))

gam = getGamma(p,j,n);

sums1 = 0;
sums2 = 0;

c = n - 2*j - 2*p;
sstar = (k - c)/2;
    
for s = 0:floor((k-1)/2)

    bcks = nchoosek(k,s);

    if odd == 0
        sums1 = sums1 - bcks*(sin((k-2*s+c)*f + c*g))*c;
        if s ~= sstar
            sums2 = sums2 + bcks*(sin((k-2*s-c)*f - c*g))*c;
        end
    else
        sums1 = sums1 + bcks*(cos((k-2*s+c)*f + c*g))*c;
        if s ~= sstar
            sums2 = sums2 + bcks*(cos((k-2*s-c)*f - c*g))*c;
        end
    end
    
end

t23 = 0;
if lk == 1
    if odd == 0
        t23 = -nchoosek(k,k/2)*sin(c*f + c*g)*c;
    else
        t23 = nchoosek(k,k/2)*cos(c*f + c*g)*c;
    end
end

sump = sump + gam*(sums1 + sums2 + t23);

end

W3fg = sump;

end


function C = getC(k,f)
C = 0;
for s = 0:floor((k-1)/2)

    C = C + 2*nchoosek(k,s)*(sin((k-2*s)*f))/(k - 2*s);
    
end
end

function Cf = getCf(k,f)
Cf = 0;
for s = 0:floor((k-1)/2)

    Cf = Cf + 2*nchoosek(k,s)*(cos((k-2*s)*f));
    
end
end

function Sf = getSf(k,f)
Sf = 0;
for s = 0:floor((k-1)/2)

    Sf = Sf - 2*nchoosek(k,s)*(sin((k-2*s)*f))*(k-2*s);
    
end
end

function KLGHt = KLGHfunc(n,j,k,a,e,eta,i,ci,G,f,g)

KLGHt = 1/G*(-2*(n + 1) - (n - 2*j)*ci/(1+ci));

end

function KLGHt = KLGHkfunc(n,j,k,a,e,eta,i,ci,G,f,g)

KLGHt = k*KLGHfunc(n,j,k,a,e,eta,i,G,ci);

end


function KLGHtn2j = KLGHn2jfunc(n,j,k,a,e,eta,i,ci,G,f,g)

KLGHtn2j = (n-2*j)*KLGHfunc(n,j,k,a,e,eta,i,G,ci);

end

function KLGHtdi2 = KLGHdi2func(n,j,k,a,e,eta,i,ci,G,f,g)

KLGHtdi2 = 1/G*((n - 2*j)*sin(i)/(1+ci)^2);

end

function KeGHt = KeGHfunc1(n,j,k,a,e,eta,i,ci,G,f,g)

KeGHt = 1/G*(-e*(2*n-1) - e*(n-2*j)*ci/(1+ci));

end

function KeGHt = KeGHefunc1(n,j,k,a,e,eta,i,ci,G,f,g)

KeGHt = 1/G*(-(2*n-1) - (n-2*j)*ci/(1+ci));

end

function KeGHt = KeGHefunc2(n,j,k,a,e,eta,i,ci,G,f,g)


KeGHt = 1/G*(-(2*n-1)*e - e*(n-2*j)*ci/(1+ci) + 2*e)*k;

end

function KeGHt = KeGHefunc3(n,j,k,a,e,eta,i,ci,G,f,g)


KeGHt = -k*(k-1)*eta^2;

end

function KeGHt = KeGHfunc1i(n,j,k,a,e,eta,i,ci,G,f,g)

KeGHt = 1/G*(-e*(2*n-1) - e*(n-2*j)*ci/(1+ci))*(n-2*j)*ci;

end

function KeGHt = KeGHfunc2(n,j,k,a,e,eta,i,ci,G,f,g)

KeGHt = 1/G*( - k*eta^2);

end

function KeGHt = KeGHgfunc2(n,j,k,a,e,eta,i,ci,G,f,g)

KeGHt = 1/G*( - k*eta^2/1);

end

function KeGHt = KeGHfunc2i(n,j,k,a,e,eta,i,ci,G,f,g)

KeGHt = 1/G*(-k*eta^2/1)*(n-2*j)*ci;

end

function KeGHt = KeGHfunc3i(n,j,k,a,e,eta,i,ci,G,f,g)

KeGHt = 1/G*(sin(i)/(1+ci)^2)*(n-2*j)*e;

end

function [K2n,K2nk1bye,K2nk2bye] = getK2nF(func, n,Jn,J2,a,e,eta,i,g,G,f,mu,Re,divSI)

ci = cos(i);

del = getDelta(n,Jn,J2,a,eta,mu,Re);

odd = 0; 

if mod(n,2) ~= 0
    % odd harmonic
    odd = 1;
end

% outer j loop
sumj1 = 0;
sumj1k1 = 0;
sumj1k2 = 0;

for j = 0:floor(n/2)
    
    if divSI == true
        beta = getBetaSP(j,n,i);
    else
        beta = getBeta(j,n,i);
    end

    % inner k loops
    sumk1 = 0;
    sumk1k1 = 0;
    sumk1k2 = 0;

    bcnj = 0;
    if odd == 0
        bcnj = nchoosek(n-2*j, n/2-j);
    end
    
    for k1 = 0:(floor(n/2)-1)
        
        if odd == 0
            % even harmonic
            k = 2*k1;
        else
            % odd harmonic
            k = 2*k1 + 1;
        end
        
        sump = 0;
        sumpk1 = 0;
        sumpk2 = 0;
        for p = 0:(floor((n-1)/2)-j)
            
            gam = getGamma(p,j,n);
            njp = n - 2*j - 2*p;
            
            if ((k - njp)/2) < 0 || ((k - njp)/2) > k
                continue;
            end
            
            bck = nchoosek(k, (k - njp)/2 );
            
            if odd == 0
                sump = sump + gam*bck*cos(njp*g);
                if k >= 1
                    sumpk1 = sumpk1 + gam*bck*cos(njp*g);
                end
                if k >= 2
                    sumpk2 = sumpk2 + gam*bck*cos(njp*g);
                end                
            else
                sump = sump + gam*bck*sin(njp*g);
                if k >=1 
                    sumpk1 = sumpk1 + gam*bck*sin(njp*g);
                end
                if k >=2 
                    sumpk2 = sumpk2 + gam*bck*sin(njp*g);
                end                
            end
        end

        alp = getAlpha(k,n,e);
        
        if odd == 0
            t12 = bcnj*nchoosek(k,k/2) + sump;
            if k >= 1
                t12k1 = bcnj*nchoosek(k,k/2) + sumpk1;
            end
            if k >= 2
                t12k2 = bcnj*nchoosek(k,k/2) + sumpk2;
            end            
        else
            t12 = sump;
            if k >= 1
                t12k1 = sumpk1;
            end
            if k >= 2
                t12k2 = sumpk2;
            end            
        end
        
        sumk1 = sumk1 + alp/2^k*(t12)*func(n,j,k,a,e,eta,i,ci,G,f,g);
        if k >= 1
            alp1 = getAlpha1(k,n,e);
            sumk1k1 = sumk1k1 + (n-k)/k*alp1/2^k*(t12k1)*func(n,j,k,a,e,eta,i,ci,G,f,g);
        end
        if k >= 2
            alp2 = getAlpha2(k,n,e);
            sumk1k2 = sumk1k2 + (n-k)/k*alp2/2^k*(t12k2)*func(n,j,k,a,e,eta,i,ci,G,f,g);
        end        
    end

    sumj1 = sumj1 + beta*sumk1;
    sumj1k1 = sumj1k1 + beta*sumk1k1;
    sumj1k2 = sumj1k2 + beta*sumk1k2;
end

K2n = del*sumj1;
K2nk1bye = del*sumj1k1;
K2nk2bye = del*sumj1k2;
end

function [dqdqcst123] = getdqdqcst123(func, n,Jn,J2,a,e,eta,i,g,G,mu,Re)

del = getDelta(n,Jn,J2,a,eta,mu,Re);

odd = 0;

if mod(n,2) ~= 0
    % odd harmonic
    odd = 1;
end

% outer j loop
sumj1 = 0;


for j = 0:floor(n/2)
    beta = getBeta(j,n,i);

    % inner k loop
    sumk1 = 0;

    bcnj = 0;
    if odd == 0
        bcnj = nchoosek(n-2*j, n/2-j);
    end
    
    for k1 = 0:(floor(n/2)-1)
        
        if odd == 0
            % even harmonic
            k = 2*k1;
        else
            % odd harmonic
            k = 2*k1 + 1;
        end

        if k < 2
            continue;
        end
        
        sump = 0;

        for p = 0:(floor((n-1)/2)-j)
            
            gam = getGamma(p,j,n);
            njp = n - 2*j - 2*p;
            
            if ((k - njp)/2) < 0 || ((k - njp)/2) > k
                continue;
            end
            
            bck = nchoosek(k, (k - njp)/2 );
            
            if odd == 0
                sump = sump + (1+eta^2*njp^2)*gam*bck*cos(njp*g);
            else
                sump = sump + (1+eta^2*njp^2)*gam*bck*sin(njp*g);
            end
        end

        alp = getAlpha2(k,n,e);
        
        if odd == 0
            t12 = bcnj*nchoosek(k,k/2) + sump;
        else
            t12 = sump;
        end
        
        spterm = -e^2*(2*n-1) - k*eta^2 -e^2*(n-2*j)*cos(i)/(1+cos(i));
        
        sumk1 = sumk1 + alp/2^k*(n-k)/k*t12*spterm;
    end

    sumj1 = sumj1 + beta*sumk1;
end

dqdqcst123 = del*sumj1;
end

function [K2nSPF,K2nSPFk1] = getK2nSPm1F(func, n,Jn,J2,a,e,eta,i,g,G,mu,Re)

del = getDelta(n,Jn,J2,a,eta,mu,Re);

odd = 0;

if mod(n,2) ~= 0
    % odd harmonic
    odd = 1;
end

% outer j loop
sumj1 = 0;
sumj1k1 = 0;

for j = 0:floor((n-2)/2)
    betasp = getBetaSPm1(j,n,i);

    % inner k loop
    sumk1 = 0;
    sumk1k1 = 0;

    bcnj = 0;
    if odd == 0
        bcnj = nchoosek(n-2*j, n/2-j);
    end
    
    for k1 = 0:(floor(n/2)-1)
        
        if odd == 0
            % even harmonic
            k = 2*k1;
        else
            % odd harmonic
            k = 2*k1 + 1;
        end
        
        sump = 0;
        sumpk1 = 0;
        for p = 0:(floor((n-1)/2)-j)
            gam = getGamma(p,j,n);
            njp = n - 2*j - 2*p;
            
            if ((k - njp)/2) < 0 || ((k - njp)/2) > k
                continue;
            end
            
            bck = nchoosek(k, (k - njp)/2 );
            
            if odd == 0
                sump = sump + gam*bck*cos(njp*g);
                if k >= 1
                    sumpk1 = sumpk1 + gam*bck*cos(njp*g);
                end
            else
                sump = sump + gam*bck*sin(njp*g);
                if k >=1 
                    sumpk1 = sumpk1 + gam*bck*sin(njp*g);
                end
            end
        end

        alp = getAlpha(k,n,e);
        
        if odd == 0
            t12 = bcnj*nchoosek(k,k/2) + sump;
            if k >= 1
                t12k1 = bcnj*nchoosek(k,k/2) + sumpk1;
            end
        else
            t12 = sump;
            if k >= 1
                t12k1 = sumpk1;
            end
        end
        sumk1 = sumk1 + alp/2^k*(t12)*func(n,j,k,a,e,eta,i,ci,G,f,g);
        if k >= 1
            alp1 = getAlpha1(k,n,e);
            sumk1k1 = sumk1k1 + (n-k)*alp1/2^k*(t12k1)*func(n,j,k,a,e,eta,i,ci,G,f,g);
        end
    end

    if isinf(betasp) == 0
        sumj1 = sumj1 + betasp*sumk1;
        sumj1k1 = sumj1k1 + betasp*sumk1k1;
    elseif sumk1 ~= 0
        disp('Inclination singularity');
    end
    
end

K2nSPF = del*sumj1;
K2nSPFk1 = del*sumj1k1;

end



function K2nSPF = getK2nSPFk1(func, n,Jn,J2,a,e,eta,i,g,G,mu,Re)

del = getDelta(n,Jn,J2,a,eta,mu,Re);

odd = 0;

if mod(n,2) ~= 0
    % odd harmonic
    odd = 1;
end

% outer j loop
sumj1 = 0;

for j = 0:floor(n/2)
    betasp = getBetaSP(j,n,i);

    % inner k loop
    sumk1 = 0;

    bcnj = 0;
    if odd == 0
        bcnj = nchoosek(n-2*j, n/2-j);
    end
    
    for k1 = 0:(floor(n/2)-1)
        
        if odd == 0
            % even harmonic
            k = 2*k1;
        else
            % odd harmonic
            k = 2*k1 + 1;
        end
        
        if k < 1
            continue;
        end
        
        sump = 0;
        for p = 0:(floor((n-1)/2)-j)
            gam = getGamma(p,j,n);
            njp = n - 2*j - 2*p;
            
            if ((k - njp)/2) < 0 || ((k - njp)/2) > k
                continue;
            end
            
            bck = nchoosek(k, (k - njp)/2 );
            
            if odd == 0
                sump = sump + gam*bck*cos(njp*g);
            else
                sump = sump + gam*bck*sin(njp*g);
            end
        end

        alp = getAlpha1(k,n,e);
        if odd == 0
            t12 = bcnj*nchoosek(k,k/2) + sump;
        else
            t12 = sump;
        end
        sumk1 = sumk1 + (n-k)*alp/2^k*(t12)*func(n,j,k,a,e,eta,i,ci,G,f,g);
    end

    if isinf(betasp) == 0
        sumj1 = sumj1 + betasp*sumk1;
    elseif sumk1 ~= 0
        disp('Inclination singularity');
    end
    
end

K2nSPF = del*sumj1;

end


function [K2ng,K2ngk1bye,K2ngk2bye] = getK2ngF(func, n,Jn,J2,a,e,eta,i,g,G,f,mu,Re,divsi)

ci = cos(i);

del = getDelta(n,Jn,J2,a,eta,mu,Re);

odd = 0;

if mod(n,2) ~= 0
    % odd harmonic
    odd = 1;
end

% outer j loop
sumj1 = 0;
sumj1k1 = 0;
sumj1k2 = 0;

for j = 0:floor(n/2)
    
    if divsi == true
        beta = getBetaSP(j,n,i);
    else
        beta = getBeta(j,n,i);
    end

    % inner k loops
    sumk1 = 0;
    sumk1k1 = 0;
    sumk1k2 = 0;

    for k1 = 0:(floor(n/2)-1)
        
        if odd == 0
            % even harmonic
            k = 2*k1;
        else
            % odd harmonic
            k = 2*k1 + 1;
        end
        
        sump = 0;
        sumpk1 = 0;
        sumpk2 = 0;
        for p = 0:(floor((n-1)/2)-j)
            
            gam = getGamma(p,j,n);
            njp = n - 2*j - 2*p;
            
            if ((k - njp)/2) < 0 || ((k - njp)/2) > k
                continue;
            end
            
            bck = nchoosek(k, (k - njp)/2 );
            
            if odd == 0
                sump = sump - gam*bck*njp*sin(njp*g);
                if k >= 1
                    sumpk1 = sumpk1 - gam*bck*njp*sin(njp*g);
                end
                if k >= 2
                    sumpk2 = sumpk2 - gam*bck*njp*sin(njp*g);
                end                
            else
                sump = sump + gam*bck*njp*cos(njp*g);
                if k >=1 
                    sumpk1 = sumpk1 + gam*bck*njp*cos(njp*g);
                end
                if k >=2 
                    sumpk2 = sumpk2 + gam*bck*njp*cos(njp*g);
                end                
            end
        end

        alp = getAlpha(k,n,e);
        
        t12 = sump;
        if k >= 1
            t12k1 = sumpk1;
        end
        if k >= 2
            t12k2 = sumpk2;
        end
        
        sumk1 = sumk1 + alp/2^k*(t12)*func(n,j,k,a,e,eta,i,ci,G,f,g);
        if k >= 1
            alp1 = getAlpha1(k,n,e);
            sumk1k1 = sumk1k1 + (n-k)/k*alp1/2^k*(t12k1)*func(n,j,k,a,e,eta,i,ci,G,f,g);
        end
        if k >= 2
            alp2 = getAlpha2(k,n,e);
            sumk1k2 = sumk1k2 + (n-k)/k*alp2/2^k*(t12k2)*func(n,j,k,a,e,eta,i,ci,G,f,g);
        end        
    end

    sumj1 = sumj1 + beta*sumk1;
    sumj1k1 = sumj1k1 + beta*sumk1k1;
    sumj1k2 = sumj1k2 + beta*sumk1k2;
end

K2ng = del*sumj1;
K2ngk1bye = del*sumj1k1;
K2ngk2bye = del*sumj1k2;

end

function [K2ng2,K2ng2k1bye,K2ng2k2bye] = getK2ng2F(func, n,Jn,J2,a,e,eta,i,g,G,f,mu,Re,divsi)

ci = cos(i);
del = getDelta(n,Jn,J2,a,eta,mu,Re);

odd = 0;

if mod(n,2) ~= 0
    % odd harmonic
    odd = 1;
end

% outer j loop
sumj1 = 0;
sumj1k1 = 0;
sumj1k2 = 0;

for j = 0:floor(n/2)
    
    if divsi == true
        beta = getBetaSP(j,n,i);
    else
        beta = getBeta(j,n,i);
    end

    % inner k loops
    sumk1 = 0;
    sumk1k1 = 0;
    sumk1k2 = 0;

    for k1 = 0:(floor(n/2)-1)
        
        if odd == 0
            % even harmonic
            k = 2*k1;
        else
            % odd harmonic
            k = 2*k1 + 1;
        end
        
        sump = 0;
        sumpk1 = 0;
        sumpk2 = 0;
        for p = 0:(floor((n-1)/2)-j)
            
            gam = getGamma(p,j,n);
            njp = n - 2*j - 2*p;
            
            if ((k - njp)/2) < 0 || ((k - njp)/2) > k
                continue;
            end
            
            bck = nchoosek(k, (k - njp)/2 );
            
            if odd == 0
                sump = sump - gam*bck*njp^2*cos(njp*g);
                if k >= 1
                    sumpk1 = sumpk1 - gam*bck*njp^2*cos(njp*g);
                end
                if k >= 2
                    sumpk2 = sumpk2 - gam*bck*njp^2*cos(njp*g);
                end                
            else
                sump = sump - gam*bck*njp^2*sin(njp*g);
                if k >=1 
                    sumpk1 = sumpk1 - gam*bck*njp^2*sin(njp*g);
                end
                if k >=2 
                    sumpk2 = sumpk2 - gam*bck*njp^2*sin(njp*g);
                end                
            end
        end

        alp = getAlpha(k,n,e);
        
        t12 = sump;
        if k >= 1
            t12k1 = sumpk1;
        end
        if k >= 2
            t12k2 = sumpk2;
        end
        
        sumk1 = sumk1 + alp/2^k*(t12)*func(n,j,k,a,e,eta,i,ci,G,f,g);
        if k >= 1
            alp1 = getAlpha1(k,n,e);
            sumk1k1 = sumk1k1 + (n-k)/k*alp1/2^k*(t12k1)*func(n,j,k,a,e,eta,i,ci,G,f,g);
        end
        if k >= 2
            alp2 = getAlpha2(k,n,e);
            sumk1k2 = sumk1k2 + (n-k)/k*alp2/2^k*(t12k2)*func(n,j,k,a,e,eta,i,ci,G,f,g);
        end        
    end

    sumj1 = sumj1 + beta*sumk1;
    sumj1k1 = sumj1k1 + beta*sumk1k1;
    sumj1k2 = sumj1k2 + beta*sumk1k2;
end

K2ng2 = del*sumj1;
K2ng2k1bye = del*sumj1k1;
K2ng2k2bye = del*sumj1k2;

end


function K2nSPg = getK2nSPgm1F(func,n,Jn,J2,a,e,eta,i,g,mu,Re)

del = getDelta(n,Jn,J2,a,eta,mu,Re);

odd = 0;

if mod(n,2) ~= 0
    % odd harmonic
    odd = 1;
end

% outer j loop
sumj1 = 0;

for j = 0:floor((n-2)/2)
    betasp = getBetaSPm1(j,n,i);

    % inner k loop
    sumk1 = 0;
    for k1 = 0:(floor(n/2)-1)
        
        if odd == 0
            % even harmonic
            k = 2*k1;
        else
            % odd harmonic
            k = 2*k1 + 1;
        end
        
        sump = 0;
        for p = 0:(floor((n-1)/2)-j)
            gam = getGamma(p,j,n);
            njp = n - 2*j - 2*p;
            
            if ((k - njp)/2) < 0 || ((k - njp)/2) > k
                continue;
            end
            
            bck = nchoosek(k, (k - njp)/2 );
            
            if odd == 0
                sump = sump - gam*bck*njp*sin(njp*g);
            else
                sump = sump + gam*bck*njp*cos(njp*g);
            end
        end

        alp = getAlpha(k,n,e);
        sumk1 = sumk1 + alp/2^k*sump*func(n,j,k,a,e,eta,i,ci,G,f,g);
    end
    
    if isinf(betasp) == 0
        sumj1 = sumj1 + betasp*sumk1;
    elseif sumk1 ~= 0
        disp('Inclination singularity');
    end
    
    
end

K2nSPg = del*sumj1;

end

function K2nSPg = getK2nSPg2m1F(func,n,Jn,J2,a,e,eta,i,g,mu,Re)

del = getDelta(n,Jn,J2,a,eta,mu,Re);

odd = 0;

if mod(n,2) ~= 0
    % odd harmonic
    odd = 1;
end

% outer j loop
sumj1 = 0;

for j = 0:floor((n-2)/2)
    betasp = getBetaSPm1(j,n,i);

    % inner k loop
    sumk1 = 0;
    for k1 = 0:(floor(n/2)-1)
        
        if odd == 0
            % even harmonic
            k = 2*k1;
        else
            % odd harmonic
            k = 2*k1 + 1;
        end
        
        sump = 0;
        for p = 0:(floor((n-1)/2)-j)
            gam = getGamma(p,j,n);
            njp = n - 2*j - 2*p;
            
            if ((k - njp)/2) < 0 || ((k - njp)/2) > k
                continue;
            end
            
            bck = nchoosek(k, (k - njp)/2 );
            
            if odd == 0
                sump = sump - gam*bck*njp^2*cos(njp*g);
            else
                sump = sump - gam*bck*njp^2*sin(njp*g);
            end
        end

        alp = getAlpha(k,n,e);
        sumk1 = sumk1 + alp/2^k*sump*func(n,j,k,a,e,eta,i,ci,G,f,g);
    end
    
    if isinf(betasp) == 0
        sumj1 = sumj1 + betasp*sumk1;
    elseif sumk1 ~= 0
        disp('Inclination singularity');
    end
    
    
end

K2nSPg = del*sumj1;

end


function K2nSPgbye = getK2nSPgbyeF(func,n,Jn,J2,a,e,eta,i,g,mu,Re)

del = getDelta(n,Jn,J2,a,eta,mu,Re);

odd = 0;

if mod(n,2) ~= 0
    % odd harmonic
    odd = 1;
end

% outer j loop
sumj1 = 0;

for j = 0:floor(n/2)
    betasp = getBetaSP(j,n,i);

    % inner k loop
    sumk1 = 0;
    for k1 = 0:(floor(n/2)-1)
        
        if odd == 0
            % even harmonic
            k = 2*k1;
        else
            % odd harmonic
            k = 2*k1 + 1;
        end
        
        if k < 1
            continue;
        end
        
        sump = 0;
        for p = 0:(floor((n-1)/2)-j)
            gam = getGamma(p,j,n);
            njp = n - 2*j - 2*p;
            
            if ((k - njp)/2) < 0 || ((k - njp)/2) > k
                continue;
            end
            
            bck = nchoosek(k, (k - njp)/2 );
            
            if odd == 0
                sump = sump - gam*bck*njp*sin(njp*g)*func(n,j,k,a,e,eta,i,ci,G,f,g);
            else
                sump = sump + gam*bck*njp*cos(njp*g)*func(n,j,k,a,e,eta,i,ci,G,f,g);
            end
        end

        alp = getAlpha1(k,n,e);
        sumk1 = sumk1 + (n-k)/k*alp/2^k*sump;
    end
    
    if isinf(betasp) == 0
        sumj1 = sumj1 + betasp*sumk1;
    elseif sumk1 ~= 0
        disp('Inclination singularity');
    end
    
    
end

K2nSPgbye = del*sumj1;

end

function [W1lge, W1lgel] = getW1lge(n,Jn,J2,a,e,eta,i,f,l,g,G,Delta,mu,Re,K2n, K2ngk1bye)

dfdl = (1 + e*cos(f))^2/eta^3;

t1 = a/(2*mu*eta)*(4*cos(f)+e*cos(2*f)+e*(2*eta^2+3*eta+3)/(1+eta))*K2n;

t1dl = a/(2*mu*eta)*(-4*sin(f)-2*e*sin(2*f))*dfdl*K2n;

t2 = -a*eta/mu*Delta*K2ngk1bye;

W1lge = t1 + t2;

W1lgel = t1dl + a*eta/mu*(dfdl-1)*K2ngk1bye;

end

function [W1lgee] = getdW1lgede(n,Jn,J2,a,e,eta,i,f,l,Delta,g,G,mu,Re,K2n,K2ne,dfde)


t1 = a/(2*mu*eta)*(4*cos(f)+e*cos(2*f)+e*(2*eta^2+3*eta+3)/(1+eta))*(e/eta^2*K2n + K2ne);

t2 = a/(2*mu*eta)*(-4*sin(f)*dfde + cos(2*f) -2*e*sin(2*f)*dfde + (4*eta^2+5*eta-1)/(1+eta))*K2n;

odd = 0;
if mod(n,2) ~= 0
    % odd harmonic
    odd = 1;
end

% outer j loop
sumj = 0; sumjk2 = 0;

for j = 0:floor(n/2)
    beta = getBeta(j,n,i);

    % inner k loop
    sumk = 0;sumk2 = 0;

    for k1 = 0:(floor(n/2)-1)
        
        if odd == 0
            % even harmonic
            k = 2*k1;
        else
            % odd harmonic
            k = 2*k1 + 1;
        end

        if k < 1
            continue;
        end
        
        sump = 0;
        for p = 0:(floor((n-1)/2)-j)
            gam = getGamma(p,j,n);
            njp = n - 2*j - 2*p;
            
            if ((k - njp)/2) < 0 || ((k - njp)/2) > k
                continue;
            end
            
            bck = nchoosek(k, (k - njp)/2 );
            
            if odd == 0
                sump = sump - gam*bck*njp*sin(njp*g);
            else
                sump = sump + gam*bck*njp*cos(njp*g);
            end
        end
        
        alp = getAlpha1(k,n,e);
        sumk = sumk + alp/2^k*(n-k)/k*sump;

        if k >= 2
            alp2 = getAlpha2(k,n,e);
            sumk2 = sumk2 + (k-1)*alp2/2^k*(n-k)/k*sump;
        end
    end
    
    sumj = sumj + beta*sumk;
    sumjk2 = sumjk2 + beta*sumk2;
    
end

del = getDelta(n,Jn,J2,a,eta,mu,Re);

t3 = a/mu*((e/eta*Delta - eta*dfde - eta*Delta*(2*n-1)*e/eta^2)*del*sumj ...
        -a*eta/mu*Delta*del*sumjk2);

W1lgee = t1 + t2 + t3;

end

function [W1lgei] = getW1lgei(n,Jn,J2,a,e,eta,i,f,l,Delta, g,G,mu,Re,K2ni)

t1si = a/(2*mu*eta)*(4*cos(f)+e*cos(2*f)+e*(2*eta^2+3*eta+3)/(1+eta))*K2ni;

odd = 0;
if mod(n,2) ~= 0
    % odd harmonic
    odd = 1;
end

% outer j loop
sumj = 0;

for j = 0:floor(n/2)
    betaSP = getBetaSP(j,n,i);

    % inner k loop
    sumk = 0;

    for k1 = 0:(floor(n/2)-1)
        
        if odd == 0
            % even harmonic
            k = 2*k1;
        else
            % odd harmonic
            k = 2*k1 + 1;
        end

        if k < 1
            continue;
        end
        
        sump = 0;
        for p = 0:(floor((n-1)/2)-j)
            gam = getGamma(p,j,n);
            njp = n - 2*j - 2*p;
            
            if ((k - njp)/2) < 0 || ((k - njp)/2) > k
                continue;
            end
            
            bck = nchoosek(k, (k - njp)/2 );
            
            if odd == 0
                sump = sump - gam*bck*njp*sin(njp*g);
            else
                sump = sump + gam*bck*njp*cos(njp*g);
            end
        end
        
        alp1 = getAlpha1(k,n,e);
        
        sumk = sumk + alp1/2^k*(n-k)/k*sump;
    end
    
    sumj = sumj + (n-2*j)*cos(i)*betaSP*sumk;
end

del = getDelta(n,Jn,J2,a,eta,mu,Re);
W1lgei = t1si - a*eta/mu*Delta*del*sumj;

end


function [W1lgegbysi] = getW1lgeg(n,Jn,J2,a,e,eta,i,f,l,Delta,g,G,mu,Re,K2ngbysi)

t1 = a/(2*mu*eta)*(4*cos(f)+e*cos(2*f)+e*(2*eta^2+3*eta+3)/(1+eta))*K2ngbysi;

odd = 0;
if mod(n,2) ~= 0
    % odd harmonic
    odd = 1;
end

% outer j loop
sumj = 0;

for j = 0:floor(n/2)
    betaSP = getBetaSP(j,n,i);

    % inner k loop
    sumk = 0;

    for k1 = 0:(floor(n/2)-1)
        
        if odd == 0
            % even harmonic
            k = 2*k1;
        else
            % odd harmonic
            k = 2*k1 + 1;
        end

        if k < 1
            continue;
        end
        
        sump = 0;
        for p = 0:(floor((n-1)/2)-j)
            gam = getGamma(p,j,n);
            njp = n - 2*j - 2*p;
            
            if ((k - njp)/2) < 0 || ((k - njp)/2) > k
                continue;
            end
            
            bck = nchoosek(k, (k - njp)/2 );
            
            if odd == 0
                sump = sump - gam*bck*njp^2*cos(njp*g);
            else
                sump = sump - gam*bck*njp^2*sin(njp*g);
            end
        end
        
        alp = getAlpha1(k,n,e);
        
        sumk = sumk + alp/2^k*(n-k)/k*sump;
    end
    
    sumj = sumj + betaSP*sumk;
    
end

del = getDelta(n,Jn,J2,a,eta,mu,Re);
W1lgegbysi = t1 - a*eta/mu*Delta*del*sumj;

end


function W23lge = getW23lge(func,n,Jn,J2,a,e,eta,i,f,l,g,L,G,mu,Re,divsi)

ci = cos(i);

odd = 0;
if mod(n,2) ~= 0
    % odd harmonic
    odd = 1;
end

dfdl = (1 + e*cos(f))^2/eta^3;

% outer j loop
sumj = 0;

for j = 0:floor(n/2)
    
    if divsi == true
        beta = getBetaSP(j,n,i);
    else
        beta = getBeta(j,n,i);
    end

    sump = 0;
    % inner p loop
    for p = 0:(floor((n-1)/2)-j)
      
        gam = getGamma(p,j,n);
        c = n - 2*j - 2*p;
       
        if odd == 0
            sump = sump + gam*cos(c*f + c*g);
        else
            sump = sump + gam*sin(c*f + c*g);
        end
    end

    t1 = 1/2/G*(4*cos(f) + e*cos(2*f) + 3*e)*sump;
    
    % inner k loop
    sumk = 0;
    for k = 1:(n-1)
        
        alp1 = getAlpha1(k,n,e);
        W3f = getW3f(n,j,k,a,e,i,eta,f,g);
        W3g = getW3g(n,j,k,a,e,i,eta,f,g);
        W3fg = W3f*dfdl - 1/eta*W3g;
        sumk = sumk + alp1/2^k*(n-k)/k*eta^2/L*W3fg;
    end

    sumj = sumj + beta*(t1 + sumk)*func(n,j,0,a,e,eta,i,ci,G,f,g);

end

del = getDeltad(n,Jn,J2,a,eta,mu,Re);
W23lge = del*sumj;

end

function W23lgee = getdW23lgede(n,Jn,J2,a,e,eta,i,f,l,g,L,G,W23lge,mu,Re,dfdl,dfde,d2fdlde)

odd = 0;
if mod(n,2) ~= 0
    % odd harmonic
    odd = 1;
end

% outer j loop
sumj = 0;

for j = 0:floor(n/2)
    beta = getBeta(j,n,i);

    sump = 0; sump2 = 0;
    
    % inner p loop
    for p = 0:(floor((n-1)/2)-j)
        gam = getGamma(p,j,n);
        c = n - 2*j - 2*p;
       
        if odd == 0
            sump = sump + gam*cos(c*f + c*g);
        else
            sump = sump + gam*sin(c*f + c*g);
        end
    end
    
    % inner p loop
    for p = 0:(floor((n-1)/2)-j)
        gam = getGamma(p,j,n);
        c = n - 2*j - 2*p;
       
        if odd == 0
            sump2 = sump2 - gam*sin(c*f + c*g)*c;
        else
            sump2 = sump2 + gam*cos(c*f + c*g)*c;
        end
    end    

    t1 = 1/2/eta/G*(e/eta)*(4*cos(f) + e*cos(2*f) + 3*e)*sump;
    t2 = 1/2/G*(-4*sin(f)*dfde + cos(2*f) - 2*e*sin(2*f)*dfde + 3)*sump;
    t3 = 1/2/G*(4*cos(f) + e*cos(2*f) + 3*e)*sump2*dfde;
    
    % inner k loop
    sumk = 0;sumk2 = 0;
    for k = 1:(n-1)
        
        alp = getAlpha1(k,n,e);
        W3f = getW3f(n,j,k,a,e,i,eta,f,g);
        W3f2 = getW3f2(n,j,k,a,e,i,eta,f,g);
        W3g = getW3g(n,j,k,a,e,i,eta,f,g);
        W3fg = getW3fg(n,j,k,a,e,i,eta,f,g);
        
        W3fgt1 = W3f*dfdl - 1/eta*W3g;
        W3fgt2 = W3f2*dfde*dfdl + W3f*d2fdlde ...
                    - 1/eta^2*e/eta*W3g -1/eta*W3fg*dfde;
                
        sumk = sumk + alp/2^k*(n-k)/k*eta^2/L*W3fgt1 ...
                + alp/2^k*(n-k)/k*eta^2/L*W3fgt2;
            
        if k >= 2
            alp2 = getAlpha2(k,n,e);
            sumk2 = sumk2 + alp2/2^k*(k-1)*(n-k)/k*eta^2/L*W3fgt1;
        end
    end
    
    sumj = sumj + beta*(t1 + t2 + t3 + sumk + sumk2);

end

del = getDeltad(n,Jn,J2,a,eta,mu,Re);
W23lgee = (2*n-1)*e/eta^2*W23lge + del*sumj;

end

function W23lgegbysi = getW23lgeg(func,n,Jn,J2,a,e,eta,i,f,l,g,L,G,dfdl,mu,Re)

ci = cos(i);

odd = 0;
if mod(n,2) ~= 0
    % odd harmonic
    odd = 1;
end

% outer j loop
sumj = 0;

for j = 0:floor(n/2)
    betaSP = getBetaSP(j,n,i);

    sump = 0;
    % inner p loop
    for p = 0:(floor((n-1)/2)-j)
        gam = getGamma(p,j,n);
        c = n - 2*j - 2*p;
       
        if odd == 0
            sump = sump - gam*sin(c*f + c*g)*c*func(n,j,0,a,e,eta,i,ci,G,f,g);
        else
            sump = sump + gam*cos(c*f + c*g)*c*func(n,j,0,a,e,eta,i,ci,G,f,g);
        end
    end

    t1 = 1/2/G*(4*cos(f) + e*cos(2*f) + 3*e)*sump;
    
    % inner k loop
    sumk = 0;
    for k = 1:(n-1)
        
        alp = getAlpha1(k,n,e);
        W3fg = getW3fg(n,j,k,a,e,i,eta,f,g);
        W3g2 = getW3g2(n,j,k,a,e,i,eta,f,g);
        W3fgterm = W3fg*dfdl - 1/eta*W3g2;
        sumk = sumk + alp/2^k*(n-k)/k*eta^2/L*W3fgterm*func(n,j,k,a,e,eta,i,ci,G,f,g);
    end
    sumj = sumj + betaSP*(t1 + sumk);

end

del = getDeltad(n,Jn,J2,a,eta,mu,Re);
W23lgegbysi = del*sumj;

end



function W23lgel = getW23lgel(n,Jn,J2,a,e,eta,i,f,l,g,L,G,mu,Re,dfdl,d2fdfdl,divsi)

odd = 0;
if mod(n,2) ~= 0
    % odd harmonic
    odd = 1;
end

% outer j loop
sumj = 0;

for j = 0:floor(n/2)
    
    if divsi == true
        beta = getBetaSP(j,n,i);
    else
        beta = getBeta(j,n,i);
    end

    sump = 0;sumpl = 0;
    % inner p loop
    for p = 0:(floor((n-1)/2)-j)
        
        gam = getGamma(p,j,n);
        c = n - 2*j - 2*p;
       
        if odd == 0
            sump = sump + gam*cos(c*f + c*g);
            sumpl = sumpl - gam*sin(c*f + c*g)*c;
        else
            sump = sump + gam*sin(c*f + c*g);
            sumpl = sumpl + gam*cos(c*f + c*g)*c;
        end
    end

    t1l = 1/2/G*(-4*sin(f) - 2*e*sin(2*f))*sump*dfdl + 1/2/G*(4*cos(f) + e*cos(2*f) + 3*e)*sumpl*dfdl;
    
    % inner k loop
    sumk = 0;
    for k = 1:(n-1)
        
        alp1 = getAlpha1(k,n,e);
        W3f = getW3f(n,j,k,a,e,i,eta,f,g);
        W3f2 = getW3f2(n,j,k,a,e,i,eta,f,g);
%         W3g = getW3g(n,j,k,a,e,i,eta,f,g);
        W3fg = getW3fg(n,j,k,a,e,i,eta,f,g);
        W3fgterm = W3f2*dfdl^2 + W3f*d2fdfdl*dfdl - 1/eta*W3fg*dfdl;
        sumk = sumk + alp1/2^k*(n-k)/k*eta^2/L*W3fgterm;
    end
    sumj = sumj + beta*(t1l + sumk);

end

del = getDeltad(n,Jn,J2,a,eta,mu,Re);
W23lgel = del*sumj;

end



function func = dummy(n,j,k,a,e,eta,i,ci,G,f,g)

func = 1;

end

function func = kterm(n,j,k,a,e,eta,i,ci,G,f,g)

func = k;

end


function func = q1q1sst2w2kterm(n,j,k,a,e,eta,i,ci,G,f,g)

func = getW3f(n,j,k,a,e,i,eta,f,g);

end
function func = q1q1sst5term1(n,j,k,a,e,eta,i,G,f,g,ci)

func = - 2*n*e - e*(n-2*j)*ci/(1+ci)  + e;

end

function func = q1q1sst5term2(n,j,k,a,e,eta,i,ci,G,f,g)

func = - 2*n - (n-2*j)*ci/(1+ci)  + (1-k*eta^2);

end

function func = q1q1cs12t1(n,j,k,a,e,eta,i,ci,G,f,g)

if mod(n,2) == 0
    ci = cos(i);
    func =  (-(2*n-1)- (n-2*j)*ci/(1+ci))*nchoosek(n-2*j,n/2-j)*nchoosek(k,k/2);
else
    func = 0;
end
end

function func = q1q1cs12t2(n,j,k,a,e,eta,i,ci,G,f,g)

ci = cos(i);

t1 =  (-e*(2*n - 1) - e*(n - 2*j)*ci/(1+ci));

% inner p loop
sump = 0;
for p = 0:(floor(((n-1)/2)-j))
    
    gam = getGamma(p,j,n);
    c = n - 2*j - 2*p;
    if odd == 0
        sump = sump + gam*cos(c*f + c*g)*func(n,j,k,a,e,eta,i,ci,G,f,g);
    else
        sump = sump + gam*sin(c*f + c*g)*func(n,j,k,a,e,eta,i,ci,G,f,g);
    end
end

func = sump*t1;

end

function func = q1q1sc5s1pcterm(n,j,k,a,e,eta,i,ci,G,f,g)

sump = 0;
for p = 0:(floor(((n-1)/2)-j))
    
    c = n - 2*j - 2*p;
    
    if mod(n,2) == 0
        tterm = sin((1+c)*f + c*g);
    else
        tterm = -cos((1+c)*f + c*g);
    end
    
    t2 =  (n-1)/2/eta^4*((-e/(1+eta)*(c*eta^4 + c*eta^3 -eta^4 -eta^3 -c*eta -eta^2 -c-eta-1) ...
        + (-2*c*eta^3 - 2*c*eta^2 + 4*c + 4)*cos(f) - e*(c*eta^3 + c*eta^2 -6*c - 6)*(cos(f))^2 ...
        + (-4*eta^2*c - 4*eta^2 + 4*c + 4)*(cos(f))^3 -e*(c*eta^2 +eta^2 -c -1)*(cos(f))^4))*tterm;
    
    sump = sump + t2;
        
end

func = sump;

end

function func = q1q1sc5s1mcterm(n,j,k,a,e,eta,i,ci,G,f,g)

sump = 0;
for p = 0:(floor(((n-1)/2)-j))
    
c = n - 2*j - 2*p;

if mod(n,2) == 0
    tterm = sin((1-c)*f - c*g);
else
    tterm = cos((1-c)*f - c*g);
end

t2 =  (n-1)/2/eta^4*((e/(1+eta)*(c*eta^4 + c*eta^3 +eta^4 +eta^3 -c*eta +eta^2 -c+eta+1) ...
            + (2*c*eta^3 + 2*c*eta^2 - 4*c + 4)*cos(f) + e*(c*eta^3 + c*eta^2 -6*c + 6)*(cos(f))^2 ...
            + (4*eta^2*c - 4*eta^2 - 4*c + 4)*(cos(f))^3 +e*(c*eta^2 -eta^2 -c +1)*(cos(f))^4))*tterm;

        sump = sump + t2;
        
end

func = sump;
        
end


function func = q1q1sc5cfgterm(n,j,k,a,e,eta,i,ci,G,f,g)

sump = 0;
for p = 0:(floor(((n-1)/2)-j))
    
    c = n - 2*j - 2*p;
    
    if mod(n,2) == 0
        tterm = cos(c*f + c*g);
    else
        tterm = sin(c*f + c*g);
    end
    sump = sump + tterm;
end

func =  sump;

end

function func = q1q1sc5cfg2term(n,j,k,a,e,eta,i,ci,G,f,g)

sump = 0;
for p = 0:(floor(((n-1)/2)-j))

c = n - 2*j - 2*p;

if mod(n,2) == 0
    tterm = -sin(c*f + c*g)*c;
else
    tterm = cos(c*f + c*g)*c;
end

sump = sump + tterm;
end

func =  sump;


end

function func = q1q1sc51pc1mcterm(n,j,k,a,e,eta,i,ci,G,f,g)
sump = 0;
for p = 0:(floor(((n-1)/2)-j))

c = n - 2*j - 2*p;

if mod(n,2) == 0
    tterm = cos((1+c)*f + c*g) + cos((1-c)*f-c*g);
else
    tterm = sin((1+c)*f + c*g) - sin((1-c)*f-c*g);
end

sump = sump + tterm;
end

func =  sump;

end

function func = q1q1sc5sfterm(n,j,k,a,e,eta,i,ci,G,f,g)

if mod(n,2) == 0
    func =  nchoosek(n-2*j,n/2-j)*(n-1)/2;
else
    func = 0;
end

end

function func = q1q1sc5cfterm(n,j,k,a,e,eta,i,ci,G,f,g)

if mod(n,2) == 0
    func =  nchoosek(n-2*j,n/2-j);
else
    func = 0;
end

end



function del = getDelta(n,Jn,J2,a,eta,mu,Re)

del = 2*Jn/J2^2*mu*Re^n/(2^n*eta^(2*n-1)*a^(n+1));

end

function deld = getDeltad(n,Jn,J2,a,eta,mu,Re)

deld = 2*Jn/J2^2*sqrt(mu)*Re^n/(2^n*eta^(2*n-1)*a^(n-1/2));

end

function beta = getBeta(j,n,i)

beta = (-1)^j*factorial(2*n-2*j)*(sin(i))^(n-2*j)/...
        (factorial(j)*factorial(n-j)*factorial(n-2*j)*2^(n-2*j));

ERR_VAL = -999999999;    
if isnan(beta) == true || isinf(beta) == true
    beta = ERR_VAL;
end
    
end

function betasp = getBetaSP(j,n,i)

betasp = (-1)^j*factorial(2*n-2*j)*(sin(i))^(n-2*j-1)/...
        (factorial(j)*factorial(n-j)*factorial(n-2*j)*2^(n-2*j));

ERR_VAL = -999999999;    
if isnan(betasp) == true || isinf(betasp) == true
    betasp = ERR_VAL;
end
    
end

function betasp = getBetaSPm1(j,n,i)

betasp = (-1)^j*factorial(2*n-2*j)*(sin(i))^(n-2*j-2)/...
        (factorial(j)*factorial(n-j)*factorial(n-2*j)*2^(n-2*j));

end

function alp = getAlpha(k,n,e)

alp = e^k*nchoosek(n-1,k);

end

function alp = getAlpha1(k,n,e)

alp = e^(k-1)*nchoosek(n-1,k-1);

end

function alp = getAlpha2(k,n,e)

alp = e^(k-2)*nchoosek(n-1,k-1);

end

function gam = getGamma(p,j,n)

s = 1;
if mod(n,2) ~= 0 
    % odd harmonic
    s = -1;
end

gam = ((-1)^(floor(n/2) - j - (s*p)))*2*nchoosek(n - 2*j, p);

end

