function [ daJ2, daJn ] = ATSecDriftCond( EqOE, dp1,dp2,dq1,dq2, mu,Re,n, Jcoeff, tol,  J22_ON, J23_ON )
%ATSecDriftCond Along-Track Secular Drift Mitigation Condition using
%Secular rates due to Zonal Harmonics, J2 up to O(3) and rest up to O(2)
% Author: Bharat Mahajan (https://github.com/princemahajan)


% Get Classical elements
RV = Eqn2RV(EqOE ,mu, tol, false);
coe = OrbitElem( RV(1:3), RV(4:6), mu );
OE = coe(2:7)';
OE(6) = coe(end);

a0 = OE(1);
e0 = OE(2);
i0 = OE(3);
h0 = OE(4);
g0 = OE(5);


de = cos(g0 + h0)*dq1 + sin(g0 + h0)*dq2;
di = (1+cos(i0))*(cos(h0)*dp1 + sin(h0)*dp2);

J2 = Jcoeff(2);

% Variations of the along-track distance y for Zonals Jn (n>2)
daJn = 0;
for ctr = 3:n
    
    Jn = Jcoeff(ctr);
    
    [ dfda, dfde, dfdi ] = ATDriftVariations(a0,e0,i0, ctr, Jn, J2, mu, Re );
    
    % condition on SMA in classical elements: dfda*da+ dfde*de + dfdi*di = 0
    if (dfda || dfde || dfdi) ~= 0
        daJn = daJn + J2^2/2*(-(dfde*de + dfdi*di)/dfda);
    end

end

% J2
if n > 1
    eta0 = sqrt(1-e0^2);
    deta = (-e0/eta0)*de;
    daJ2 = ATJ2DRiftVariations(a0,e0,i0,deta,di, J2, mu, Re, J22_ON, J23_ON );
else
    daJ2 = 0;
end

% % condition on SMA in classical elements: dfda*da+ dfde*de + dfdi*di = 0
% da = daJn + daJ2;


end



function [ dfda, dfde, dfdi ] = ATDriftVariations(a,e,i, n, Jn, J2, mu, Re )

dfda = 0;
dfde = 0;
dfdi = 0;

% No contributions from odd zonal harmonics
if mod(n,2) ~= 0
    return;
end

% intermediate terms
eta = sqrt(1 - e^2);
L = sqrt(mu*a);
% G = L*eta;
% H = G*cos(i);
si = sin(i);
ci = cos(i);
% ti = si/ci;

dLda = 1/2*sqrt(mu/a);

sumhdotsi = 0;
sumdfda = 0;
sumdfde = 0;
sumdfdi = 0;

for j = 0:floor(n/2)
    
    bn2 = nchoosek((n-2*j),(n/2-j));
    
    sumhdotsik = 0;
    sumdfdak = 0;
    sumdfdek = 0;
    
    for k1 = 0:(floor(n/2)-1)
        
        if mod(n,2) == 0
            % even harmonic
            k = 2*k1;
        else
            % odd harmonic
            k = 2*k1 + 1;
        end
        
        Ckk2 = nchoosek(k,k/2);
        [al, ale] = alpha(k,n,e);
        
        t1 = al/2^k*Ckk2;
        t1e = ale/2^k*Ckk2;

        % dfda
        sumdfdak = sumdfdak + t1*(-3 -(2*n - 1)/eta - k*eta/(1 + eta));

        % dfde
        sumdfdek = sumdfdek + t1e*(-3 -(2*n - 1)/eta - k*eta/(1 + eta)) ...
             + t1*((2*n - 1)/eta^2*(-e/eta) - k*((1+eta)*(-e/eta)-eta*(-e/eta))/(1 + eta)^2);

        % hdot
        if (n-2*j) ~= 0
            sumhdotsik = sumhdotsik + t1*(-1/eta*(n-2*j)*ci/si);
        else
            sumhdotsik = sumhdotsik + 0;
        end
    end
    
    [bt,bti] = betajn(j,n,i);

    sumdfda = sumdfda + bt*bn2*sumdfdak;
    sumdfde = sumdfde + bt*bn2*sumdfdek;
    sumdfdi = sumdfdi + bti*bn2*sumdfdak;
    
    sumhdotsi = sumhdotsi + bt*bn2*sumhdotsik;
    
end

[del, dela, dele] = delta(n, a, eta, mu, Re, J2, Jn);

% Variation w.r.t. a
dfda = (1/L*dela + del*(-1/L^2)*dLda)*sumdfda;

% Variation w.r.t. e
dfde = 1/L*dele*sumdfda + 1/L*del*sumdfde;

% Variation w.r.t. i
dfdi = 1/L*del*sumdfdi + 1/L*del*sumhdotsi;

end


function [del, dela, dele] = delta(n, a, eta, mu, Re, J2, Jn)

del = 2*Jn/J2^2*mu*Re^n/(2^n*eta^(2*n-1)*a^(n+1));

% derivatove w.r.t. a
dela = 2*Jn/J2^2*mu*Re^n/(2^n*eta^(2*n-1))*(-n-1)*a^(-n-2);

% derivatove w.r.t. e
e = sqrt(1-eta^2);
dele = 2*Jn/J2^2*mu*Re^n/(2^n*a^(n+1))*(-2*n+1)*eta^(-2*n)*(-e/eta);

end


function [b, bi] = betajn(j,n,i)

n2j = n - 2*j;

b = (-1)^j*factorial(2*n-2*j)*((sin(i))^(n2j))/(factorial(j)*factorial(n-j)*factorial(n2j)*2^(n2j));

% derivative w.r.t. i
bi = (-1)^j*factorial(2*n-2*j)/(factorial(j)*factorial(n-j)*factorial(n2j)*2^(n2j))*n2j*((sin(i))^(n2j-1))*cos(i);

if n2j == 0
    bi = 0;
end

end



function [al,ale] = alpha(k,n,e)

al = e^k*(nchoosek((n-1), (k)));

% derivative w.r.t. e
if k == 0
    ale = 0;
else
    ale = k*e^(k-1)*(nchoosek((n-1), (k)));
end

end

function [da] = ATJ2DRiftVariations(a,e,i,deta,di, J2, mu, Re, J22_ON, J23_ON )

J = J2;

R__e = Re;

eta = sqrt(1 - e^2);

da = -0.2e1 / 0.3e1 * J * R__e ^ 2 * a * ((0.3e1 * R__e ^ 4 * (0.7172e4 - 0.111650000e9 * cos(i) ^ 12 + (147 * eta ^ 10) + 0.245388e6 * cos(i) ^ 2 - 0.5671160e7 * cos(i) ^ 4 + 0.40525320e8 * cos(i) ^ 6 - 0.129277500e9 * cos(i) ^ 8 + 0.188699500e9 * cos(i) ^ 10 + (94365 * eta ^ 2) - (258414 * eta ^ 4) - (148796 * eta ^ 5) + (183928 * eta ^ 6) + (52173 * eta) + 0.543320000e9 * cos(i) ^ 6 * (eta ^ 4) + 0.653649120e9 * cos(i) ^ 6 * (eta ^ 3) - 0.13681834e8 * cos(i) ^ 4 * (eta ^ 5) + 0.463727490e9 * cos(i) ^ 6 * (eta ^ 2) - 0.96641190e8 * cos(i) ^ 4 * (eta ^ 4) + 0.203201130e9 * cos(i) ^ 6 * eta - 0.99987970e8 * cos(i) ^ 4 * (eta ^ 3) + 0.2306110e7 * cos(i) ^ 2 * (eta ^ 5) - 0.63471090e8 * cos(i) ^ 4 * (eta ^ 2) + 0.8245920e7 * cos(i) ^ 2 * (eta ^ 4) - 0.27277470e8 * cos(i) ^ 4 * eta + 0.6084432e7 * cos(i) ^ 2 * (eta ^ 3) + 0.2565315e7 * cos(i) ^ 2 * (eta ^ 2) + 0.982047e6 * cos(i) ^ 2 * eta + 0.41861460e8 * cos(i) ^ 6 * (eta ^ 5) - 0.3656625e7 * cos(i) ^ 10 * (eta ^ 10) + 0.2772175e7 * cos(i) ^ 8 * (eta ^ 10) - 0.889350e6 * cos(i) ^ 6 * (eta ^ 10) + 0.127400e6 * cos(i) ^ 4 * (eta ^ 10) - 0.7497e4 * cos(i) ^ 2 * (eta ^ 10) + 0.32838750e8 * cos(i) ^ 12 * (eta ^ 9) + 0.204911250e9 * cos(i) ^ 12 * (eta ^ 8) + 0.1653750e7 * cos(i) ^ 12 * (eta ^ 10) + 0.564363750e9 * cos(i) ^ 12 * (eta ^ 7) + 0.609016250e9 * cos(i) ^ 12 * (eta ^ 6) - 0.151431250e9 * cos(i) ^ 12 * (eta ^ 5) - 0.1246341250e10 * cos(i) ^ 12 * (eta ^ 4) - 0.1764433750e10 * cos(i) ^ 12 * (eta ^ 3) - 0.1377030000e10 * cos(i) ^ 12 * (eta ^ 2) - 0.598537500e9 * cos(i) ^ 12 * eta + (10183 * eta ^ 9) + (84290 * eta ^ 8) + (241810 * eta ^ 7) - (40298 * eta ^ 3) + 0.78588290e8 * cos(i) ^ 4 * (eta ^ 7) - 0.326373e6 * cos(i) ^ 2 * (eta ^ 9) - 0.1507484475e10 * cos(i) ^ 8 * (eta ^ 2) + 0.78815994e8 * cos(i) ^ 4 * (eta ^ 6) - 0.2493060e7 * cos(i) ^ 2 * (eta ^ 8) - 0.663097875e9 * cos(i) ^ 8 * eta - 0.7150280e7 * cos(i) ^ 2 * (eta ^ 7) - 0.6581282e7 * cos(i) ^ 2 * (eta ^ 6) - 0.73822125e8 * cos(i) ^ 10 * (eta ^ 9) - 0.464158500e9 * cos(i) ^ 10 * (eta ^ 8) - 0.1277557000e10 * cos(i) ^ 10 * (eta ^ 7) + 0.60041875e8 * cos(i) ^ 8 * (eta ^ 9) - 0.1378507250e10 * cos(i) ^ 10 * (eta ^ 6) + 0.387113150e9 * cos(i) ^ 8 * (eta ^ 8) + 0.167248750e9 * cos(i) ^ 10 * (eta ^ 5) + 0.1073158150e10 * cos(i) ^ 8 * (eta ^ 7) - 0.21918830e8 * cos(i) ^ 6 * (eta ^ 9) + 0.2161796000e10 * cos(i) ^ 10 * (eta ^ 4) + 0.1149804500e10 * cos(i) ^ 8 * (eta ^ 6) - 0.147175000e9 * cos(i) ^ 6 * (eta ^ 8) + 0.2949314000e10 * cos(i) ^ 10 * (eta ^ 3) - 0.89920200e8 * cos(i) ^ 8 * (eta ^ 5) - 0.412454960e9 * cos(i) ^ 6 * (eta ^ 7) + 0.3897416e7 * cos(i) ^ 4 * (eta ^ 9) + 0.2255314875e10 * cos(i) ^ 10 * (eta ^ 2) - 0.1552307050e10 * cos(i) ^ 8 * (eta ^ 4) - 0.432231660e9 * cos(i) ^ 6 * (eta ^ 6) + 0.27693934e8 * cos(i) ^ 4 * (eta ^ 8) + 0.988032375e9 * cos(i) ^ 10 * eta - 0.2023199550e10 * cos(i) ^ 8 * (eta ^ 3)) * deta + 0.4e1 * cos(i) * sin(i) * R__e ^ 4 * eta * (0.48132e5 - (1134 * eta ^ 10) - 0.971400e6 * cos(i) ^ 2 + 0.7854840e7 * cos(i) ^ 4 - 0.31240800e8 * cos(i) ^ 6 + 0.60598500e8 * cos(i) ^ 8 - 0.45675000e8 * cos(i) ^ 10 + (551041 * eta ^ 2) + (602903 * eta ^ 4) - (140741 * eta ^ 5) - (777247 * eta ^ 6) + (239427 * eta) - 0.370922650e9 * cos(i) ^ 6 * (eta ^ 4) - 0.508277150e9 * cos(i) ^ 6 * (eta ^ 3) - 0.16219470e8 * cos(i) ^ 4 * (eta ^ 5) - 0.372806000e9 * cos(i) ^ 6 * (eta ^ 2) + 0.96807450e8 * cos(i) ^ 4 * (eta ^ 4) - 0.161212500e9 * cos(i) ^ 6 * eta + 0.128707710e9 * cos(i) ^ 4 * (eta ^ 3) + 0.2327763e7 * cos(i) ^ 2 * (eta ^ 5) + 0.92295510e8 * cos(i) ^ 4 * (eta ^ 2) - 0.12304695e8 * cos(i) ^ 2 * (eta ^ 4) + 0.39925170e8 * cos(i) ^ 4 * eta - 0.15998685e8 * cos(i) ^ 2 * (eta ^ 3) - 0.11268666e8 * cos(i) ^ 2 * (eta ^ 2) - 0.4876704e7 * cos(i) ^ 2 * eta + 0.52505050e8 * cos(i) ^ 6 * (eta ^ 5) + 0.1063125e7 * cos(i) ^ 10 * (eta ^ 10) - 0.1779750e7 * cos(i) ^ 8 * (eta ^ 10) + 0.1021650e7 * cos(i) ^ 6 * (eta ^ 10) - 0.270900e6 * cos(i) ^ 4 * (eta ^ 10) + 0.32865e5 * cos(i) ^ 2 * (eta ^ 10) - (33426 * eta ^ 9) - (227855 * eta ^ 8) - (655001 * eta ^ 7) + (785581 * eta ^ 3) - 0.94156230e8 * cos(i) ^ 4 * (eta ^ 7) + 0.712041e6 * cos(i) ^ 2 * (eta ^ 9) + 0.737786625e9 * cos(i) ^ 8 * (eta ^ 2) - 0.109951530e9 * cos(i) ^ 4 * (eta ^ 6) + 0.4535499e7 * cos(i) ^ 2 * (eta ^ 8) + 0.318463875e9 * cos(i) ^ 8 * eta + 0.12697665e8 * cos(i) ^ 2 * (eta ^ 7) + 0.14838477e8 * cos(i) ^ 2 * (eta ^ 6) + 0.19288125e8 * cos(i) ^ 10 * (eta ^ 9) + 0.115381875e9 * cos(i) ^ 10 * (eta ^ 8) + 0.316475625e9 * cos(i) ^ 10 * (eta ^ 7) - 0.32744250e8 * cos(i) ^ 8 * (eta ^ 9) + 0.371698125e9 * cos(i) ^ 10 * (eta ^ 6) - 0.197178375e9 * cos(i) ^ 8 * (eta ^ 8) + 0.13411875e8 * cos(i) ^ 10 * (eta ^ 5) - 0.541046625e9 * cos(i) ^ 8 * (eta ^ 7) + 0.19608450e8 * cos(i) ^ 6 * (eta ^ 9) - 0.510819375e9 * cos(i) ^ 10 * (eta ^ 4) - 0.634986375e9 * cos(i) ^ 8 * (eta ^ 6) + 0.120033250e9 * cos(i) ^ 6 * (eta ^ 8) - 0.743173125e9 * cos(i) ^ 10 * (eta ^ 3) - 0.68152125e8 * cos(i) ^ 8 * (eta ^ 5) + 0.331148950e9 * cos(i) ^ 6 * (eta ^ 7) - 0.5449500e7 * cos(i) ^ 4 * (eta ^ 9) - 0.571008750e9 * cos(i) ^ 10 * (eta ^ 2) + 0.694215375e9 * cos(i) ^ 8 * (eta ^ 4) + 0.387687350e9 * cos(i) ^ 6 * (eta ^ 6) - 0.33927690e8 * cos(i) ^ 4 * (eta ^ 8) - 0.245362500e9 * cos(i) ^ 10 * eta + 0.982135125e9 * cos(i) ^ 8 * (eta ^ 3)) * di) * J22_ON*J ^ 2 + (-0.48e2 * R__e ^ 2 * a ^ 2 * (eta ^ 4) * (0.280e3 - 0.4760e4 * cos(i) ^ 2 + 0.27440e5 * cos(i) ^ 4 - 0.47600e5 * cos(i) ^ 6 - 0.77000e5 * cos(i) ^ 8 + 0.245000e6 * cos(i) ^ 10 + (2239 * eta ^ 2) - (2315 * eta ^ 4) - (3745 * eta ^ 5) - (2543 * eta ^ 6) + (1337 * eta) + 0.2329150e7 * cos(i) ^ 6 * (eta ^ 4) + 0.1237850e7 * cos(i) ^ 6 * (eta ^ 3) - 0.591810e6 * cos(i) ^ 4 * (eta ^ 5) + 0.160810e6 * cos(i) ^ 6 * (eta ^ 2) - 0.495510e6 * cos(i) ^ 4 * (eta ^ 4) - 0.136570e6 * cos(i) ^ 6 * eta - 0.112290e6 * cos(i) ^ 4 * (eta ^ 3) + 0.75265e5 * cos(i) ^ 2 * (eta ^ 5) + 0.143486e6 * cos(i) ^ 4 * (eta ^ 2) + 0.53315e5 * cos(i) ^ 2 * (eta ^ 4) + 0.118258e6 * cos(i) ^ 4 * eta - 0.3815e4 * cos(i) ^ 2 * (eta ^ 3) - 0.34079e5 * cos(i) ^ 2 * (eta ^ 2) - 0.22057e5 * cos(i) ^ 2 * eta + 0.2268650e7 * cos(i) ^ 6 * (eta ^ 5) - (125 * eta ^ 8) - (871 * eta ^ 7) + (815 * eta ^ 3) - 0.117454e6 * cos(i) ^ 4 * (eta ^ 7) - 0.2326925e7 * cos(i) ^ 8 * (eta ^ 2) - 0.362462e6 * cos(i) ^ 4 * (eta ^ 6) + 0.2325e4 * cos(i) ^ 2 * (eta ^ 8) - 0.653275e6 * cos(i) ^ 8 * eta + 0.16431e5 * cos(i) ^ 2 * (eta ^ 7) + 0.49023e5 * cos(i) ^ 2 * (eta ^ 6) + 0.15625e5 * cos(i) ^ 10 * (eta ^ 8) + 0.204875e6 * cos(i) ^ 10 * (eta ^ 7) + 0.1070875e7 * cos(i) ^ 10 * (eta ^ 6) - 0.65625e5 * cos(i) ^ 8 * (eta ^ 8) + 0.3073125e7 * cos(i) ^ 10 * (eta ^ 5) - 0.543675e6 * cos(i) ^ 8 * (eta ^ 7) + 0.5379375e7 * cos(i) ^ 10 * (eta ^ 4) - 0.2002275e7 * cos(i) ^ 8 * (eta ^ 6) + 0.51250e5 * cos(i) ^ 6 * (eta ^ 8) + 0.5908125e7 * cos(i) ^ 10 * (eta ^ 3) - 0.4230125e7 * cos(i) ^ 8 * (eta ^ 5) + 0.385910e6 * cos(i) ^ 6 * (eta ^ 7) + 0.3981125e7 * cos(i) ^ 10 * (eta ^ 2) - 0.5551375e7 * cos(i) ^ 8 * (eta ^ 4) + 0.1262230e7 * cos(i) ^ 6 * (eta ^ 6) - 0.16250e5 * cos(i) ^ 4 * (eta ^ 8) + 0.1505875e7 * cos(i) ^ 10 * eta - 0.4596125e7 * cos(i) ^ 8 * (eta ^ 3)) * deta - 0.192e3 * cos(i) * sin(i) * R__e ^ 2 * a ^ 2 * (eta ^ 5) * (-0.40e2 + 0.320e3 * cos(i) ^ 2 + 0.1200e4 * cos(i) ^ 4 - 0.16000e5 * cos(i) ^ 6 + 0.35000e5 * cos(i) ^ 8 + (47 * eta ^ 2) + (1735 * eta ^ 4) + (1775 * eta ^ 5) + (1029 * eta ^ 6) - (131 * eta) - 0.729500e6 * cos(i) ^ 6 * (eta ^ 4) - 0.650500e6 * cos(i) ^ 6 * (eta ^ 3) + 0.193950e6 * cos(i) ^ 4 * (eta ^ 5) - 0.362800e6 * cos(i) ^ 6 * (eta ^ 2) + 0.232650e6 * cos(i) ^ 4 * (eta ^ 4) - 0.115400e6 * cos(i) ^ 6 * eta + 0.172350e6 * cos(i) ^ 4 * (eta ^ 3) - 0.30680e5 * cos(i) ^ 2 * (eta ^ 5) + 0.74910e5 * cos(i) ^ 4 * (eta ^ 2) - 0.32860e5 * cos(i) ^ 2 * (eta ^ 4) + 0.16530e5 * cos(i) ^ 4 * eta - 0.19940e5 * cos(i) ^ 2 * (eta ^ 3) - 0.5464e4 * cos(i) ^ 2 * (eta ^ 2) + 0.208e3 * cos(i) ^ 2 * eta - 0.526000e6 * cos(i) ^ 6 * (eta ^ 5) + (45 * eta ^ 8) + (327 * eta ^ 7) + (845 * eta ^ 3) + 0.29010e5 * cos(i) ^ 4 * (eta ^ 7) + 0.594875e6 * cos(i) ^ 8 * (eta ^ 2) + 0.99330e5 * cos(i) ^ 4 * (eta ^ 6) - 0.700e3 * cos(i) ^ 2 * (eta ^ 8) + 0.219625e6 * cos(i) ^ 8 * eta - 0.5204e4 * cos(i) ^ 2 * (eta ^ 7) - 0.16912e5 * cos(i) ^ 2 * (eta ^ 6) + 0.3125e4 * cos(i) ^ 8 * (eta ^ 8) + 0.37375e5 * cos(i) ^ 8 * (eta ^ 7) + 0.184625e6 * cos(i) ^ 8 * (eta ^ 6) - 0.7500e4 * cos(i) ^ 6 * (eta ^ 8) + 0.506875e6 * cos(i) ^ 8 * (eta ^ 5) - 0.63300e5 * cos(i) ^ 6 * (eta ^ 7) + 0.854375e6 * cos(i) ^ 8 * (eta ^ 4) - 0.239400e6 * cos(i) ^ 6 * (eta ^ 6) + 0.3750e4 * cos(i) ^ 4 * (eta ^ 8) + 0.908125e6 * cos(i) ^ 8 * (eta ^ 3)) * di) * J - 0.1536e4 * a ^ 4 * (eta ^ 8) * (0.4e1 - 0.72e2 * cos(i) ^ 2 + 0.480e3 * cos(i) ^ 4 - 0.1400e4 * cos(i) ^ 6 + 0.1500e4 * cos(i) ^ 8 + (55 * eta ^ 2) + (50 * eta ^ 4) + (19 * eta ^ 5) + (3 * eta ^ 6) + (23 * eta) - 0.17500e5 * cos(i) ^ 6 * (eta ^ 4) - 0.24500e5 * cos(i) ^ 6 * (eta ^ 3) + 0.2280e4 * cos(i) ^ 4 * (eta ^ 5) - 0.19250e5 * cos(i) ^ 6 * (eta ^ 2) + 0.6000e4 * cos(i) ^ 4 * (eta ^ 4) - 0.8050e4 * cos(i) ^ 6 * eta + 0.8400e4 * cos(i) ^ 4 * (eta ^ 3) - 0.342e3 * cos(i) ^ 2 * (eta ^ 5) + 0.6600e4 * cos(i) ^ 4 * (eta ^ 2) - 0.900e3 * cos(i) ^ 2 * (eta ^ 4) + 0.2760e4 * cos(i) ^ 4 * eta - 0.1260e4 * cos(i) ^ 2 * (eta ^ 3) - 0.990e3 * cos(i) ^ 2 * (eta ^ 2) - 0.414e3 * cos(i) ^ 2 * eta - 0.6650e4 * cos(i) ^ 6 * (eta ^ 5) + (70 * eta ^ 3) + 0.20625e5 * cos(i) ^ 8 * (eta ^ 2) + 0.360e3 * cos(i) ^ 4 * (eta ^ 6) + 0.8625e4 * cos(i) ^ 8 * eta - 0.54e2 * cos(i) ^ 2 * (eta ^ 6) + 0.1125e4 * cos(i) ^ 8 * (eta ^ 6) + 0.7125e4 * cos(i) ^ 8 * (eta ^ 5) + 0.18750e5 * cos(i) ^ 8 * (eta ^ 4) - 0.1050e4 * cos(i) ^ 6 * (eta ^ 6) + 0.26250e5 * cos(i) ^ 8 * (eta ^ 3)) * deta - 0.3072e4 * cos(i) * sin(i) * a ^ 4 * (eta ^ 9) * (-0.4e1 + 0.60e2 * cos(i) ^ 2 - 0.300e3 * cos(i) ^ 4 + 0.500e3 * cos(i) ^ 6 - (55 * eta ^ 2) - (50 * eta ^ 4) - (19 * eta ^ 5) - (3 * eta ^ 6) - (23 * eta) + 0.6250e4 * cos(i) ^ 6 * (eta ^ 4) + 0.8750e4 * cos(i) ^ 6 * (eta ^ 3) - 0.1425e4 * cos(i) ^ 4 * (eta ^ 5) + 0.6875e4 * cos(i) ^ 6 * (eta ^ 2) - 0.3750e4 * cos(i) ^ 4 * (eta ^ 4) + 0.2875e4 * cos(i) ^ 6 * eta - 0.5250e4 * cos(i) ^ 4 * (eta ^ 3) + 0.285e3 * cos(i) ^ 2 * (eta ^ 5) - 0.4125e4 * cos(i) ^ 4 * (eta ^ 2) + 0.750e3 * cos(i) ^ 2 * (eta ^ 4) - 0.1725e4 * cos(i) ^ 4 * eta + 0.1050e4 * cos(i) ^ 2 * (eta ^ 3) + 0.825e3 * cos(i) ^ 2 * (eta ^ 2) + 0.345e3 * cos(i) ^ 2 * eta + 0.2375e4 * cos(i) ^ 6 * (eta ^ 5) - (70 * eta ^ 3) - 0.225e3 * cos(i) ^ 4 * (eta ^ 6) + 0.45e2 * cos(i) ^ 2 * (eta ^ 6) + 0.375e3 * cos(i) ^ 6 * (eta ^ 6)) * di) / eta / (0.5e1 * R__e ^ 6 * (0.1793e4 - 0.27912500e8 * cos(i) ^ 12 + (63 * eta ^ 10) + 0.61347e5 * cos(i) ^ 2 - 0.1417790e7 * cos(i) ^ 4 + 0.10131330e8 * cos(i) ^ 6 - 0.32319375e8 * cos(i) ^ 8 + 0.47174875e8 * cos(i) ^ 10 + (22499 * eta ^ 2) - (77934 * eta ^ 4) - (29784 * eta ^ 5) + (77134 * eta ^ 6) + (13414 * eta) + 0.143494020e9 * cos(i) ^ 6 * (eta ^ 4) + 0.174496800e9 * cos(i) ^ 6 * (eta ^ 3) + 0.1153844e7 * cos(i) ^ 4 * (eta ^ 5) + 0.118777242e9 * cos(i) ^ 6 * (eta ^ 2) - 0.26080460e8 * cos(i) ^ 4 * (eta ^ 4) + 0.50813340e8 * cos(i) ^ 6 * eta - 0.27231700e8 * cos(i) ^ 4 * (eta ^ 3) + 0.227860e6 * cos(i) ^ 2 * (eta ^ 5) - 0.16352792e8 * cos(i) ^ 4 * (eta ^ 2) + 0.2312390e7 * cos(i) ^ 2 * (eta ^ 4) - 0.6794860e7 * cos(i) ^ 4 * eta + 0.1749424e7 * cos(i) ^ 2 * (eta ^ 3) + 0.680295e6 * cos(i) ^ 2 * (eta ^ 2) + 0.239946e6 * cos(i) ^ 2 * eta - 0.15696840e8 * cos(i) ^ 6 * (eta ^ 5) - 0.1567125e7 * cos(i) ^ 10 * (eta ^ 10) + 0.1188075e7 * cos(i) ^ 8 * (eta ^ 10) - 0.381150e6 * cos(i) ^ 6 * (eta ^ 10) + 0.54600e5 * cos(i) ^ 4 * (eta ^ 10) - 0.3213e4 * cos(i) ^ 2 * (eta ^ 10) + 0.12757500e8 * cos(i) ^ 12 * (eta ^ 9) + 0.74997500e8 * cos(i) ^ 12 * (eta ^ 8) + 0.708750e6 * cos(i) ^ 12 * (eta ^ 10) + 0.200783500e9 * cos(i) ^ 12 * (eta ^ 7) + 0.228747500e9 * cos(i) ^ 12 * (eta ^ 6) + 0.22500e5 * cos(i) ^ 12 * (eta ^ 5) - 0.322575000e9 * cos(i) ^ 12 * (eta ^ 4) - 0.460837500e9 * cos(i) ^ 12 * (eta ^ 3) - 0.351790250e9 * cos(i) ^ 12 * (eta ^ 2) - 0.150550000e9 * cos(i) ^ 12 * eta + (3858 * eta ^ 9) + (30205 * eta ^ 8) + (85748 * eta ^ 7) - (19476 * eta ^ 3) + 0.27871908e8 * cos(i) ^ 4 * (eta ^ 7) - 0.124398e6 * cos(i) ^ 2 * (eta ^ 9) - 0.385327905e9 * cos(i) ^ 8 * (eta ^ 2) + 0.30259392e8 * cos(i) ^ 4 * (eta ^ 6) - 0.898345e6 * cos(i) ^ 2 * (eta ^ 8) - 0.166154250e9 * cos(i) ^ 8 * eta - 0.2533664e7 * cos(i) ^ 2 * (eta ^ 7) - 0.2593306e7 * cos(i) ^ 2 * (eta ^ 6) - 0.28662750e8 * cos(i) ^ 10 * (eta ^ 9) - 0.169772625e9 * cos(i) ^ 10 * (eta ^ 8) - 0.454605600e9 * cos(i) ^ 10 * (eta ^ 7) + 0.23258250e8 * cos(i) ^ 8 * (eta ^ 9) - 0.517742250e9 * cos(i) ^ 10 * (eta ^ 6) + 0.141298925e9 * cos(i) ^ 8 * (eta ^ 8) - 0.43943500e8 * cos(i) ^ 10 * (eta ^ 5) + 0.381561820e9 * cos(i) ^ 8 * (eta ^ 7) - 0.8457780e7 * cos(i) ^ 6 * (eta ^ 9) + 0.561133750e9 * cos(i) ^ 10 * (eta ^ 4) + 0.432825350e9 * cos(i) ^ 8 * (eta ^ 6) - 0.53545350e8 * cos(i) ^ 6 * (eta ^ 8) + 0.774590000e9 * cos(i) ^ 10 * (eta ^ 3) + 0.48558400e8 * cos(i) ^ 8 * (eta ^ 5) - 0.146470848e9 * cos(i) ^ 6 * (eta ^ 7) + 0.1495656e7 * cos(i) ^ 4 * (eta ^ 9) + 0.576171775e9 * cos(i) ^ 10 * (eta ^ 2) - 0.405122350e9 * cos(i) ^ 8 * (eta ^ 4) - 0.163791420e9 * cos(i) ^ 6 * (eta ^ 6) + 0.10031898e8 * cos(i) ^ 4 * (eta ^ 8) + 0.248020250e9 * cos(i) ^ 10 * eta - 0.534681500e9 * cos(i) ^ 8 * (eta ^ 3)) * J23_ON*J ^ 3 - 0.176e3 * R__e ^ 4 * a ^ 2 * (eta ^ 4) * (0.35e2 - 0.595e3 * cos(i) ^ 2 + 0.3430e4 * cos(i) ^ 4 - 0.5950e4 * cos(i) ^ 6 - 0.9625e4 * cos(i) ^ 8 + 0.30625e5 * cos(i) ^ 10 + (264 * eta ^ 2) - (450 * eta ^ 4) - (670 * eta ^ 5) - (464 * eta ^ 6) + (166 * eta) + 0.382500e6 * cos(i) ^ 6 * (eta ^ 4) + 0.203700e6 * cos(i) ^ 6 * (eta ^ 3) - 0.102940e6 * cos(i) ^ 4 * (eta ^ 5) + 0.34560e5 * cos(i) ^ 6 * (eta ^ 2) - 0.86500e5 * cos(i) ^ 4 * (eta ^ 4) - 0.15260e5 * cos(i) ^ 6 * eta - 0.26180e5 * cos(i) ^ 4 * (eta ^ 3) + 0.13310e5 * cos(i) ^ 2 * (eta ^ 5) + 0.14736e5 * cos(i) ^ 4 * (eta ^ 2) + 0.9850e4 * cos(i) ^ 2 * (eta ^ 4) + 0.14444e5 * cos(i) ^ 4 * eta + 0.1010e4 * cos(i) ^ 2 * (eta ^ 3) - 0.3904e4 * cos(i) ^ 2 * (eta ^ 2) - 0.2726e4 * cos(i) ^ 2 * eta + 0.385100e6 * cos(i) ^ 6 * (eta ^ 5) - (25 * eta ^ 8) - (166 * eta ^ 7) + (30 * eta ^ 3) - 0.22284e5 * cos(i) ^ 4 * (eta ^ 7) - 0.323800e6 * cos(i) ^ 8 * (eta ^ 2) - 0.65376e5 * cos(i) ^ 4 * (eta ^ 6) + 0.465e3 * cos(i) ^ 2 * (eta ^ 8) - 0.86450e5 * cos(i) ^ 8 * eta + 0.3126e4 * cos(i) ^ 2 * (eta ^ 7) + 0.8904e4 * cos(i) ^ 2 * (eta ^ 6) + 0.3125e4 * cos(i) ^ 10 * (eta ^ 8) + 0.36750e5 * cos(i) ^ 10 * (eta ^ 7) + 0.177000e6 * cos(i) ^ 10 * (eta ^ 6) - 0.13125e5 * cos(i) ^ 8 * (eta ^ 8) + 0.473750e6 * cos(i) ^ 10 * (eta ^ 5) - 0.101550e6 * cos(i) ^ 8 * (eta ^ 7) + 0.781250e6 * cos(i) ^ 10 * (eta ^ 4) - 0.349200e6 * cos(i) ^ 8 * (eta ^ 6) + 0.10250e5 * cos(i) ^ 6 * (eta ^ 8) + 0.816250e6 * cos(i) ^ 10 * (eta ^ 3) - 0.691750e6 * cos(i) ^ 8 * (eta ^ 5) + 0.72860e5 * cos(i) ^ 6 * (eta ^ 7) + 0.528000e6 * cos(i) ^ 10 * (eta ^ 2) - 0.856250e6 * cos(i) ^ 8 * (eta ^ 4) + 0.225040e6 * cos(i) ^ 6 * (eta ^ 6) - 0.3250e4 * cos(i) ^ 4 * (eta ^ 8) + 0.193250e6 * cos(i) ^ 10 * eta - 0.672250e6 * cos(i) ^ 8 * (eta ^ 3)) * J22_ON*J ^ 2 - 0.3584e4 * R__e ^ 2 * a ^ 4 * (eta ^ 8) * (0.1e1 - 0.18e2 * cos(i) ^ 2 + 0.120e3 * cos(i) ^ 4 - 0.350e3 * cos(i) ^ 6 + 0.375e3 * cos(i) ^ 8 + (15 * eta ^ 2) + (15 * eta ^ 4) + (6 * eta ^ 5) + (eta ^ 6) + (6 * eta) - 0.5250e4 * cos(i) ^ 6 * (eta ^ 4) - 0.7000e4 * cos(i) ^ 6 * (eta ^ 3) + 0.720e3 * cos(i) ^ 4 * (eta ^ 5) - 0.5250e4 * cos(i) ^ 6 * (eta ^ 2) + 0.1800e4 * cos(i) ^ 4 * (eta ^ 4) - 0.2100e4 * cos(i) ^ 6 * eta + 0.2400e4 * cos(i) ^ 4 * (eta ^ 3) - 0.108e3 * cos(i) ^ 2 * (eta ^ 5) + 0.1800e4 * cos(i) ^ 4 * (eta ^ 2) - 0.270e3 * cos(i) ^ 2 * (eta ^ 4) + 0.720e3 * cos(i) ^ 4 * eta - 0.360e3 * cos(i) ^ 2 * (eta ^ 3) - 0.270e3 * cos(i) ^ 2 * (eta ^ 2) - 0.108e3 * cos(i) ^ 2 * eta - 0.2100e4 * cos(i) ^ 6 * (eta ^ 5) + (20 * eta ^ 3) + 0.5625e4 * cos(i) ^ 8 * (eta ^ 2) + 0.120e3 * cos(i) ^ 4 * (eta ^ 6) + 0.2250e4 * cos(i) ^ 8 * eta - 0.18e2 * cos(i) ^ 2 * (eta ^ 6) + 0.375e3 * cos(i) ^ 8 * (eta ^ 6) + 0.2250e4 * cos(i) ^ 8 * (eta ^ 5) + 0.5625e4 * cos(i) ^ 8 * (eta ^ 4) - 0.350e3 * cos(i) ^ 6 * (eta ^ 6) + 0.7500e4 * cos(i) ^ 8 * (eta ^ 3)) * J - 0.2048e4 * a ^ 6 * (eta ^ 12) * (0.125e3 * cos(i) ^ 6 * (eta ^ 5) + 0.625e3 * cos(i) ^ 6 * (eta ^ 4) + 0.1250e4 * cos(i) ^ 6 * (eta ^ 3) - 0.75e2 * cos(i) ^ 4 * (eta ^ 5) + 0.1250e4 * cos(i) ^ 6 * (eta ^ 2) - 0.375e3 * cos(i) ^ 4 * (eta ^ 4) + 0.625e3 * cos(i) ^ 6 * eta - 0.750e3 * cos(i) ^ 4 * (eta ^ 3) + 0.15e2 * cos(i) ^ 2 * (eta ^ 5) + 0.125e3 * cos(i) ^ 6 - 0.750e3 * cos(i) ^ 4 * (eta ^ 2) + 0.75e2 * cos(i) ^ 2 * (eta ^ 4) - 0.375e3 * cos(i) ^ 4 * eta + 0.150e3 * cos(i) ^ 2 * (eta ^ 3) - (eta ^ 5) - 0.75e2 * cos(i) ^ 4 + 0.150e3 * cos(i) ^ 2 * (eta ^ 2) - (5 * eta ^ 4) + 0.75e2 * cos(i) ^ 2 * eta - (10 * eta ^ 3) + 0.15e2 * cos(i) ^ 2 - (10 * eta ^ 2) - (5 * eta) - 0.1e1));



end