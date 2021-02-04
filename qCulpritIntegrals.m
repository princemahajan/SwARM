function [I, dIdX, Idldg1, Idldg2,d2I1dX2,d2I2dX2,dI1lgdX,dI2lgdX] = qCulpritIntegrals(coeI,f0,n,m,p,delta,ddeltada,d2deltada2,tol, HessianOn)
%qCulpritIntegrals Culprit Integrals and its partial derivatives using
%Numerical quadrature. Uses mean anoamaly as the independent variable for
%computing Jacobians and Hessians. Quadrature is evaluated w.r.t. true
%anomaly however.
% Author: Bharat Mahajan (https://github.com/princemahajan)

% classical elements

% a = coeI(1);
e = coeI(2);
TA = coeI(3);

% values at f = 0 using eccentricity expansions
[ I0,dIdX0, Idldg10, Idldg20,I1xx0,I2xx0,I1lgx0,I2lgx0 ] = eCulpritIntegrals( delta,ddeltada,d2deltada2,e,0,0,n,m,p, HessianOn );

eta = sqrt(1-e^2);

if TA < 0
    period = -2*pi;
else
    period = 2*pi;
end

E = @(f) (mod(2*atan(sqrt((1-e)/(1+e))*tan(f/2)),period) + fix(f/period)*period);
M = @(f) E(f) - e*sin(E(f));

dfde = @(f) (2+e*cos(f)).*sin(f)/eta^2;
dfdl = @(f) (1 + e*cos(f))^2/eta^3;
% dldf = @(f) 1/dfdl(f);

% Integrand argument
A = @(f) (n-2*p)*f - m*delta*M(f);

dAda = @(f) -m*ddeltada*M(f); 
dAde = @(f) (n-2*p)*dfde(f);
dAdf = (n-2*p);
dAdl = -m*delta;

% Integrands 
I1_term = @(f) cos(A(f)).*((1+e*cos(f)).^(n-1));
I2_term = @(f) sin(A(f)).*((1+e*cos(f)).^(n-1));

%%%%%%%%%%%%% Culprit Integrals %%%%%%%%%%%%%
I1 = integral(@(f) I1_term(f), f0, TA, 'RelTol',tol,'AbsTol',tol*1e-3);
I2 = integral(@(f) I2_term(f), f0, TA, 'RelTol',tol,'AbsTol',tol*1e-3);

% I1 and I2
I = [I1; I2] + I0;

%%%%%%%%%%%%% Culprit Integrals Jacobian and Hessian %%%%%%%%%%%%%

% Default values
d2I1dX2 = zeros(6)*NaN;
d2I2dX2 = zeros(6)*NaN;
dI1lgdX = zeros(2,6)*NaN;
dI2lgdX = zeros(2,6)*NaN;

%%%%%%%%%%%%% Partials w.r.t. M %%%%%%%%%%%%%

dI1df = I1_term(TA);
dI2df = I2_term(TA);

dIdl = [dI1df; dI2df]*dfdl(TA);

if HessianOn == true
    
    % d/da(d/dl I)
    
    I1la = -I2_term(TA)*dAda(TA)*dfdl(TA);
    I2la = I1_term(TA)*dAda(TA)*dfdl(TA);
    
    % d/de(d/dl I)
    
%     I1le = 3*e/eta^2*I1_term(TA)*dfdl(TA) - I2_term(TA)*dAde(TA)*dfdl(TA) + cos(A(TA))*(n+1)*(1+e*cos(TA))^n/eta^3*(cos(TA)-e*sin(TA).*dfde(TA));
%     I2le = 3*e/eta^2*I2_term(TA)*dfdl(TA) + I1_term(TA)*dAde(TA)*dfdl(TA) + sin(A(TA))*(n+1)*(1+e*cos(TA))^n/eta^3*(cos(TA)-e*sin(TA).*dfde(TA));
    
    I1le = 3*e/eta^2*I1_term(TA)*dfdl(TA) - I2_term(TA)*dAde(TA)*dfdl(TA) + (n+1)*I1_term(TA)*dfdl(TA)/eta^2*(e*cos(TA)^2+cos(TA)-2*e);
    I2le = 3*e/eta^2*I2_term(TA)*dfdl(TA) + I1_term(TA)*dAde(TA)*dfdl(TA) + (n+1)*I2_term(TA)*dfdl(TA)/eta^2*(e*cos(TA)^2+cos(TA)-2*e);
    
    
    % d/dl(d/dl I)
    
    I1ll = -I2_term(TA)*dfdl(TA)*(dAdl + dAdf*dfdl(TA)) + dfdl(TA)*cos(A(TA))*(n+1)*((1+e*cos(TA))^n)/eta^3*(-e*sin(TA));
    I2ll = I1_term(TA)*dfdl(TA)*(dAdl + dAdf*dfdl(TA))  + dfdl(TA)*sin(A(TA))*(n+1)*((1+e*cos(TA))^n)/eta^3*(-e*sin(TA));
    
end

%%%%%%%%%%%%% Partials w.r.t. SMA %%%%%%%%%%%%%

dI1da = m*ddeltada*integral(@(f) M(f).*I2_term(f),f0,TA,'RelTol',tol,'AbsTol',tol*1e-3);
dI2da = -m*ddeltada*integral(@(f) M(f).*I1_term(f),f0,TA,'RelTol',tol,'AbsTol',tol*1e-3);

dIda = [dI1da; dI2da] + dIdX0(:,1);

% (d/de of inetgrands of I1*dfdl and I2*dfdl)/(1+e*cos(f))^(n-1)
dI1intdfdlde = @(f) 1/eta^2*(-3/2*e*(n-1)*cos(A(f)) + e/4*(2*p+1)*cos(A(f)-2*f) + (-n/2+2*p+1/2)*cos(A(f)-f) ...
                                + (3/2*n-2*p+1/2)*cos(A(f)+f) + e/4*(2*n-2*p+1)*cos(A(f)+2*f));

dI2intdfdlde = @(f) 1/eta^2*(-3/2*e*(n-1)*sin(A(f)) + e/4*(2*p+1)*sin(A(f)-2*f) + (-n/2+2*p+1/2)*sin(A(f)-f) ...
                                + (3/2*n-2*p+1/2)*sin(A(f)+f) + e/4*(2*n-2*p+1)*sin(A(f)+2*f));

if HessianOn == true
    
    % d2Idada
    
    I1aa = I1xx0(1,1) + d2deltada2/ddeltada*dI1da ...
        + m*ddeltada*integral(@(f) M(f).*I1_term(f).*dAda(f),f0,TA,'RelTol',tol,'AbsTol',tol*1e-3);
    
    I2aa = I2xx0(1,1) + d2deltada2/ddeltada*dI2da ...
        - m*ddeltada*integral(@(f) -M(f).*I2_term(f).*dAda(f),f0,TA,'RelTol',tol,'AbsTol',tol*1e-3);
    
    % d2Idade
    
    I1ae = I1xx0(1,2) + m*ddeltada*integral(@(f) M(f).*dI2intdfdlde(f).*(1+e*cos(f)).^(n-1),f0,TA,'RelTol',tol,'AbsTol',tol*1e-3);
    
    I2ae = I2xx0(1,2) - m*ddeltada*integral(@(f) M(f).*dI1intdfdlde(f).*(1+e*cos(f)).^(n-1),f0,TA,'RelTol',tol,'AbsTol',tol*1e-3);
    
    % d2Idadl
    
    I1al = m*ddeltada*M(TA)*I2_term(TA)*dfdl(TA);
    I2al = -m*ddeltada*M(TA)*I1_term(TA)*dfdl(TA);
    
end

%%%%%%%%%%%%% Partials w.r.t. e %%%%%%%%%%%%%

dI1de = integral(@(f) dI1intdfdlde(f).*(1+e*cos(f)).^(n-1), f0, TA, 'RelTol',tol,'AbsTol',tol*1e-3);

dI2de = integral(@(f) dI2intdfdlde(f).*(1+e*cos(f)).^(n-1), f0, TA, 'RelTol',tol,'AbsTol',tol*1e-3);

dIde = [dI1de; dI2de] + dIdX0(:,2);

if HessianOn == true
    
    % d2Ideda
    
    I1eabar = integral(@(f) -dI2intdfdlde(f).*dAda(f).*(1+e*cos(f)).^(n-1), f0, TA, 'RelTol',tol,'AbsTol',tol*1e-3);

    I1ea = I1xx0(2,1) + I1eabar;
    
    I2eabar = integral(@(f) dI1intdfdlde(f).*dAda(f).*(1+e*cos(f)).^(n-1), f0, TA, 'RelTol',tol,'AbsTol',tol*1e-3);

    I2ea = I2xx0(2,1) + I2eabar;
    
    % d2Idede
    
    Ieeterm1 = @(f) (-3/2*1*(n-1)*cos(A(f)) - (-3/2*e*(n-1)*sin(A(f)).*dAde(f)) ...
                    + 1/4*(2*p+1)*cos(A(f)-2*f) - e/4*(2*p+1)*sin(A(f)-2*f).*(dAde(f)-2*dfde(f)) ...
                    - (-n/2+2*p+1/2)*sin(A(f)-f).*(dAde(f)-dfde(f))  ...
                    - (3/2*n-2*p+1/2)*sin(A(f)+f).*(dAde(f)+dfde(f))  ...
                    + 1/4*(2*n-2*p+1)*cos(A(f)+2*f)) - e/4*(2*n-2*p+1)*sin(A(f)+2*f).*(dAde(f)+2*dfde(f)) ...
                    + dI1intdfdlde(f)*(n+1).*(e*cos(f).^2 + cos(f)-2*e);

    Ieeterm2 = @(f) (-3/2*1*(n-1)*sin(A(f)) + (-3/2*e*(n-1)*cos(A(f)).*dAde(f)) ...
                    + 1/4*(2*p+1)*sin(A(f)-2*f) + e/4*(2*p+1)*cos(A(f)-2*f).*(dAde(f)-2*dfde(f)) ...
                    + (-n/2+2*p+1/2)*cos(A(f)-f).*(dAde(f)-dfde(f))  ...
                    + (3/2*n-2*p+1/2)*cos(A(f)+f).*(dAde(f)+dfde(f))  ...
                    + 1/4*(2*n-2*p+1)*sin(A(f)+2*f)) + e/4*(2*n-2*p+1)*cos(A(f)+2*f).*(dAde(f)+2*dfde(f)) ...
                    + dI2intdfdlde(f)*(n+1).*(e*cos(f).^2 + cos(f)-2*e);
                
    I1eebar = integral(@(f) Ieeterm1(f).*(1+e*cos(f)).^(n-1), f0, TA, 'RelTol',tol,'AbsTol',tol*1e-3);

    I1ee = I1xx0(2,2) + 5*e/eta^2*dI1de + 1/eta^2*I1eebar;
    
    
    I2eebar = integral(@(f) Ieeterm2(f).*(1+e*cos(f)).^(n-1), f0, TA, 'RelTol',tol,'AbsTol',tol*1e-3);

    I2ee = I2xx0(2,2) + 5*e/eta^2*dI2de + 1/eta^2*I2eebar;
    
    % d2Idedl
    
    I1el = dI1intdfdlde(TA)*dfdl(TA)*(1+e*cos(TA))^(n-1);
    
    I2el = dI2intdfdlde(TA)*dfdl(TA)*(1+e*cos(TA))^(n-1);   

    
end

%%%%%%%%%%%%% Classical elements Results %%%%%%%%%%%%%

dIdX = [dIda, dIde, zeros(2,3), dIdl];

if HessianOn == true
    
    d2I1dX2 = [I1aa,I1ae,0,0,0,I1al;
        I1ea,I1ee,0,0,0,I1el;
        0,0,0,0,0,0;
        0,0,0,0,0,0;
        0,0,0,0,0,0;
        I1la,I1le,0,0,0,I1ll];
    
    d2I2dX2 = [I2aa,I2ae,0,0,0,I2al;
        I2ea,I2ee,0,0,0,I2el;
        0,0,0,0,0,0;
        0,0,0,0,0,0;
        0,0,0,0,0,0;
        I2la,I2le,0,0,0,I2ll];
    
end

%%%%%%%%%%%%% Partials of terms for Equinoctial Elements %%%%%%%%%%%%%

%%%%%%%%%%%%% 1/e*(I1*(dfdl-1/eta)) %%%%%%%%%%%%%

dfdlm1byetabye = @(f) 1/eta^3*(e + 2*cos(f) + e*(cos(f)).^2);

I1dldg1 = integral(@(f) I1_term(f).*dfdlm1byetabye(f), f0, TA, 'RelTol',tol,'AbsTol',tol*1e-3);
I2dldg1 = integral(@(f) I2_term(f).*dfdlm1byetabye(f), f0, TA, 'RelTol',tol,'AbsTol',tol*1e-3);

Idldg1 = [I1dldg1; I2dldg1] + Idldg10;

if HessianOn == true
    
    % Partials w.r.t. a
    
    I1lg1a = I1lgx0(1,1) + integral(@(f) -I2_term(f).*dAda(f).*dfdlm1byetabye(f), f0, TA, 'RelTol',tol,'AbsTol',tol*1e-3);
    I2lg1a = I2lgx0(1,1) + integral(@(f) I1_term(f).*dAda(f).*dfdlm1byetabye(f), f0, TA, 'RelTol',tol,'AbsTol',tol*1e-3);
    
    % Partials w.r.t. e
    
    dfdlm1byetabyee = @(f) -3/eta^4*(-e/eta)*(e + 2*cos(f) + e*(cos(f)).^2) + 1/eta^3*(1 + (cos(f)).^2) ...
                            + 1/eta^3*(-2*sin(f).*dfde(f) - 2*e*cos(f).*sin(f).*dfde(f));
    
    I1lg1e = I1lgx0(1,2) + integral(@(f) dI1intdfdlde(f).*dfdlm1byetabye(f).*(1+e*cos(f)).^(n-1) ...
                    + I1_term(f).*dfdlm1byetabyee(f), f0, TA, 'RelTol',tol,'AbsTol',tol*1e-3);
    
    I2lg1e = I2lgx0(1,2) + integral(@(f) dI2intdfdlde(f).*dfdlm1byetabye(f).*(1+e*cos(f)).^(n-1) ...
                    + I2_term(f).*dfdlm1byetabyee(f), f0, TA, 'RelTol',tol,'AbsTol',tol*1e-3);
    
    % Partials w.r.t. l
    
    I1lg1l = I1_term(TA).*dfdl(TA).*dfdlm1byetabye(TA);
    I2lg1l = I2_term(TA).*dfdl(TA).*dfdlm1byetabye(TA);
    
    dI1lg1dX = [I1lg1a,I1lg1e,0,0,0,I1lg1l];
    
    dI2lg1dX = [I2lg1a,I2lg1e,0,0,0,I2lg1l];
    
end


%%%%%%%%%%%%% I1*sin(f)*(1+e*cos(f)) %%%%%%%%%%%%%

I1dldg2 = integral(@(f) I1_term(f).*sin(f).*(1+e*cos(f)), f0, TA, 'RelTol',tol,'AbsTol',tol*1e-3);
I2dldg2 = integral(@(f) I2_term(f).*sin(f).*(1+e*cos(f)), f0, TA, 'RelTol',tol,'AbsTol',tol*1e-3);

Idldg2 = [I1dldg2; I2dldg2] + Idldg20;

if HessianOn == true
    
    % Partials w.r.t. a
    
    I1lg2a = I1lgx0(2,1) + integral(@(f) -I2_term(f).*dAda(f).*sin(f).*(1+e*cos(f)), f0, TA, 'RelTol',tol,'AbsTol',tol*1e-3);
    I2lg2a = I2lgx0(2,1) + integral(@(f) I1_term(f).*dAda(f).*sin(f).*(1+e*cos(f)), f0, TA, 'RelTol',tol,'AbsTol',tol*1e-3);
    
    % Partials w.r.t. e
    
    dsf1ecfde = @(f) sin(f).*cos(f) + cos(f).*(1+e*cos(f)).*dfde(f) + sin(f).*(-e*sin(f)).*dfde(f);
    
    I1lg2e = I1lgx0(2,2) + integral(@(f) dI1intdfdlde(f).*sin(f).*(1+e*cos(f)).^n ...
                                + I1_term(f).*dsf1ecfde(f) , f0, TA, 'RelTol',tol,'AbsTol',tol*1e-3);
    
    I2lg2e = I2lgx0(2,2) + integral(@(f) dI2intdfdlde(f).*sin(f).*(1+e*cos(f)).^n ...
                                + I2_term(f).*dsf1ecfde(f) , f0, TA, 'RelTol',tol,'AbsTol',tol*1e-3);
    
    % Partials w.r.t. l
    
    I1lg2l = I1_term(TA).*sin(TA).*(1+e*cos(TA)).*dfdl(TA);
    I2lg2l = I2_term(TA).*sin(TA).*(1+e*cos(TA)).*dfdl(TA);
    
    dI1lg2dX = [I1lg2a,I1lg2e,0,0,0,I1lg2l];
    dI2lg2dX = [I2lg2a,I2lg2e,0,0,0,I2lg2l];
    
    dI1lgdX = [dI1lg1dX; dI1lg2dX];
    dI2lgdX = [dI2lg1dX; dI2lg2dX];
    
end

end

