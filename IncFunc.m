%IncFunc Inclination Function using the Allan form. 
% See Yan, H., Vadali, S. R., & Alfriend, K. T. (2014). A Simplified Formulation
% of the Satellite Perturbed Relative Motion Problem.
%
% Author: Bharat Mahajan (https://github.com/princemahajan)


function [ F, Fd, Fbysi2, Fdd, Fbysi2d] = IncFunc( n,m,p,i )

% coefficient free of summation
t1 = (-1)^fix((n-m+1)/2)*factorial(n+m)/((2^n)*factorial(n-p)*factorial(p));

si2 = sin(i/2);
ci2 = cos(i/2);

% j range
ji = max([0, n-m-2*p]);
jf = min([n-m, 2*n-2*p]);

% j summation
jsum = 0;
jsumbysi2 = 0;
jsumd = 0;
jsumd2 = 0;
jsumbysi2d = 0;

for jctr = ji:1:jf
    
    % binomial coeff-1
    try
        n1 = 2*n-2*p;
        k1 = jctr;
        bc1 = nchoosek(n1,k1);
    catch
        bc1 = 0;
    end
    
    % binomial coeff-2
    try
        n2 = 2*p;
        k2 = n-m-jctr;
        bc2 = nchoosek(n2,k2);
    catch
        bc2 = 0;
    end    

    % sum for F
    jsum = jsum + ((-1)^jctr)*bc1*bc2*(si2.^(m-n+2*p+2*jctr)).*ci2.^(3*n-m-2*jctr-2*p);
    % sum for F/s(i/2)
    jsumbysi2 = jsumbysi2 + ((-1)^jctr)*bc1*bc2*(si2.^(m-n+2*p+2*jctr-1)).*ci2.^(3*n-m-2*jctr-2*p);
    
    % First Derivative
    if (m-n+2*p+2*jctr) == 0
        sterm = 0;
    else
        sterm = (m-n+2*p+2*jctr)*(si2.^(m-n+2*p+2*jctr-1)).*ci2/2.*(ci2.^(3*n-m-2*jctr-2*p));
    end
    if (3*n-m-2*jctr-2*p) == 0
        cterm = 0;
    else
        cterm = -(si2.^(m-n+2*p+2*jctr))*(3*n-m-2*jctr-2*p).*(ci2.^(3*n-m-2*jctr-2*p-1)).*si2/2;
    end
    
    jsumd = jsumd + ((-1)^jctr)*bc1*bc2*(sterm + cterm);
    
    % derivative of F/si2
    if (m-n+2*p+2*jctr-1) == 0
        sterm = 0;
    else
        sterm = (m-n+2*p+2*jctr-1)*(si2.^(m-n+2*p+2*jctr-2)).*ci2/2.*(ci2.^(3*n-m-2*jctr-2*p));
    end
    if (3*n-m-2*jctr-2*p) == 0
        cterm = 0;
    else
        cterm = -(si2.^(m-n+2*p+2*jctr-1))*(3*n-m-2*jctr-2*p).*(ci2.^(3*n-m-2*jctr-2*p-1)).*si2/2;
    end

    jsumbysi2d = jsumbysi2d + ((-1)^jctr)*bc1*bc2*(sterm + cterm);
    
    % Second Derivative
    s = m - n + 2*p + 2*jctr;
    
    csterm1 = 0;
    csterm2 = 0;
    if s ~= 0 && s ~= 1
        csterm2 = s*(s-1)*ci2.^(2*n-s+2)*si2.^(s-2);
    end
    if (2*n-s) ~= 0 && (2*n-s) ~= 1
        csterm1 = (2*n-s)*(2*n-s-1)*si2.^(s+2).*ci2.^(2*n-s-2);
    end
    
    jsumd2 = jsumd2 + ((-1)^jctr)*bc1*bc2*1/4*(ci2.^(2*n-s).*si2.^s.*(-(2*n-s)*(s+1) - (2*n-s+1)*s) + csterm1 + csterm2);
    
end

Fbysi2 = t1*jsumbysi2;
F = t1*jsum;

% Derivatives of Inclination Function
Fd = t1*jsumd;
Fdd = t1*jsumd2;
Fbysi2d = t1*jsumbysi2d;

end

