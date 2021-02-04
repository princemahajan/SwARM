function [W, dWdX, didGWgdidHWh, dWdletadWdg, d2WdX2,Wghx,Wlgx] = TesseralQGFr(coe,f0,n,m,C20,Cnm,Snm,mu,we,Re,tol,HessianOn,QuadTesseralsOn)
% Closed-form Short-period and m-daily tesseral generating function in the
% Earth-fixed frame. Uses mean anoamaly as the independent variable for
% computing Jacobians and Hessians.
% Author: Bharat Mahajan (https://github.com/princemahajan)

% classical elements
a = coe(1);
e = coe(2);
i = coe(3);
hr = coe(4);
g = coe(5);
f = coe(6);

eta = sqrt(1 - e^2);
L = sqrt(mu*a);
G = L*eta;

mm = sqrt(mu/a^3);
detade = -e/eta;

delta = we/mm;
ddeltada = we/sqrt(mu)*3/2*sqrt(a);
d2deltada2 = we/sqrt(mu)*3/4/sqrt(a);

% dfdl = (1 + e*cos(f))^2/eta^3;
% dldf = 1/dfdl;
% dfde = (2+e*cos(f)).*sin(f)/eta^2;

% d2ldf2 = 2*e*eta^3*sin(f)/(1+e*cos(f))^3;
% d2ldedf = ((1+e*cos(f))^2*3*eta^2*detade - eta^3*2*(1+e*cos(f))*(cos(f)-e*sin(f)*dfde))/(1+e*cos(f))^4;
% 
% d2fdlde = -3/eta^4*(-e/eta)*(1+e*cos(f))^2 + 2*(1+e*cos(f))/eta^3*(cos(f)-e*sin(f)*dfde);
% d2fdldf = 2*(1+e*cos(f))/eta^3*(-e*sin(f));

% select proper C and S order
if mod(n-m,2) == 0
    CS1 = Cnm;
    CS2 = Snm;
else
    CS1 = -Snm;
    CS2 = Cnm;
end

% intermediate terms

C = -2/C20^2*Re^n*sqrt(mu)/a^(n-1/2)/eta^(2*n-1);
dCda = -2/C20^2*Re^n*sqrt(mu)*(-n+1/2)*a^(-n+1/2-1)/eta^(2*n-1);
dCde = -2/C20^2*Re^n*sqrt(mu)/a^(n-1/2)*(-2*n+1)*eta^(-2*n)*detade;

Caa = -2/C20^2*Re^n*sqrt(mu)*(-n+1/2)*(-n+1/2-1)*a^(-n+1/2-2)/eta^(2*n-1);
Cae = -2/C20^2*Re^n*sqrt(mu)*(-n+1/2)*a^(-n+1/2-1)*(-2*n+1)/eta^(2*n)*detade;

Cea = Cae;
Cee = -2/C20^2*Re^n*sqrt(mu)/a^(n-1/2)*(-2*n+1)*(-eta^(-2*n-3)*(2*n*e^2 + 1));

dCdX = [dCda, dCde, zeros(1,4)];
d2CdX2 = [Caa, Cae, 0,0,0,0;
          Cea, Cee, 0,0,0,0;
          zeros(4,6)];

W = 0;
dWdX = zeros(1,6);
d2WdX2 = zeros(6);

didGWgdidHWh = 0;
dWdletadWdg = 0;

Wghx = zeros(1,6);
Wlgx = zeros(1,6);

% P Summation
for p = 0:n
    
    if f < 0
        period = -2*pi;
    else
        period = 2*pi;
    end
    
    E = mod(2*atan(sqrt((1-e)/(1+e))*tan(f/2)),period) + fix(f/period)*period;
    l = E - e*sin(E);

    % beta and partials
    beta = (n-2*p)*g + m*(hr + delta*l);
    cbeta = cos(beta);
    sbeta = sin(beta);

    dbetada = m*ddeltada*l;
    dbetadg = n-2*p;
    dbetadh = m;
    dbetadl = m*delta;

    betax = [dbetada, 0, 0, dbetadh,dbetadg, dbetadl];
    betaxx = [m*d2deltada2*l, 0, 0, 0, 0, m*ddeltada;
                zeros(1,6);
                zeros(1,6);
                zeros(1,6);
                zeros(1,6);
                m*ddeltada, 0,0,0,0,0];
    
    % The culprit integral and partials
    coeI = [a,e,f];

    if QuadTesseralsOn == true
        [I, dIdX, dIdlg1, dIdlg2,I1xx,I2xx,dI1lgdX,dI2lgdX] = qCulpritIntegrals(coeI,f0,n,m,p,delta,ddeltada,d2deltada2,tol, HessianOn);
    else
        [I, dIdX, dIdlg1, dIdlg2,I1xx,I2xx,dI1lgdX,dI2lgdX] = eCulpritIntegrals( delta,ddeltada,d2deltada2,e,f,l,n,m,p,HessianOn);
    end
    
    I1x = dIdX(1,:);
    I2x = dIdX(2,:);

    % Inclination function
    [ F, Fd, Fbysi2, Fdd, Fbysi2d] = IncFunc( n,m,p,i );
    dFdX = [0,0,Fd,0,0,0];
    d2FdX2 = [zeros(2,6);
              0,0,Fdd,0,0,0;
              zeros(3,6)];    
    % W
    Vp = (CS1*(cbeta*I(1) - sbeta*(I(2))) + CS2*(cbeta*I(2) + sbeta*I(1)));
    W = W + C*F*Vp;

    % Partials of W w.r.t classical elements
    Vpx = @(Y) CS1*(cbeta*(I1x(Y)-betax(Y)*I(2)) - sbeta*(I2x(Y) + betax(Y)*I(1))) ...
        + CS2*(cbeta*(betax(Y)*I(1) + I2x(Y)) + sbeta*(I1x(Y) - betax(Y)*I(2)));

    for ctr = 1:6
        dWdX(ctr) = dWdX(ctr) + dCdX(ctr)*F*Vp + C*dFdX(ctr)*Vp + C*F*Vpx(ctr); 
    end
    
    % didG*dWdg+didH*dWdh term
    if (n-2*p) == m
        Wghcoeff = -m*tan(i/2)/G;
        Fform = F;
    else
        Wghcoeff = 1/(G*2*cos(i/2))*(cos(i)*(n-2*p)-m);
        Fform = Fbysi2;
    end
    Vgh = CS1*(-cbeta*I(2) - sbeta*(I(1))) + CS2*(cbeta*I(1) - sbeta*I(2));
    didGWgdidHWh = didGWgdidHWh + C*Wghcoeff*Fform*Vgh;
    
    % 1/e*(dWdl - 1/eta*dWdg) term
    Vlg = (CS1*(cbeta*(-(n-2*p)*dIdlg1(2) - (n+1)/eta^3*dIdlg2(1)) ...
        + sbeta*(-(n-2*p)*dIdlg1(1) + (n+1)/eta^3*dIdlg2(2))) ...
        + CS2*(cbeta*((n-2*p)*dIdlg1(1) - (n+1)/eta^3*dIdlg2(2)) ...
        + sbeta*(-(n-2*p)*dIdlg1(2) - (n+1)/eta^3*dIdlg2(1))));
    dWdletadWdg = dWdletadWdg + C*F*Vlg;
    
    
    if HessianOn == true

        dI1lgadX = dI1lgdX(1,:);
        dI1lgbdX = dI1lgdX(2,:);
        dI2lgadX = dI2lgdX(1,:);
        dI2lgbdX = dI2lgdX(2,:);
        
        Vpxy = @(X,Y) CS1*(cbeta*(I1xx(X,Y) - betax(X)*(betax(Y)*I(1) + I2x(Y)) - betax(Y)*I2x(X) - betaxx(X,Y)*I(2)) ...
                           - sbeta*(I2xx(X,Y) + betax(X)*(-betax(Y)*I(2) + I1x(Y)) + betax(Y)*I1x(X) + betaxx(X,Y)*I(1))) ...
                    + CS2*(cbeta*(I2xx(X,Y) + betax(Y)*(-betax(X)*I(2) + I1x(X)) + betax(X)*I1x(Y) + betaxx(X,Y)*I(1)) ...
                           + sbeta*(I1xx(X,Y) - betax(Y)*(betax(X)*I(1) + I2x(X)) - betax(X)*I2x(Y) - betaxx(X,Y)*I(2)));
        
        FVpx = @(X) dFdX(X)*Vp + F*Vpx(X);
        FVpxy = @(X,Y) d2FdX2(X,Y)*Vp + dFdX(X)*Vpx(Y) + dFdX(Y)*Vpx(X) + F*Vpxy(X,Y);
        
        % compute each element of Hessian of W
        for row = 1:6
            for col = 1:6
                d2WdX2(row,col) = d2WdX2(row,col) + d2CdX2(row,col)*F*Vp + dCdX(row)*FVpx(col) ...
                                            + dCdX(col)*FVpx(row) + C*FVpxy(row,col);
            end
        end
        
        % Partials of terms for equinoctial variations
        % didG*dWdg+didH*dWdh term
        if (n-2*p) == m
            Wghcoeffx = [m*tan(i/2)/(2*a*G), -m*tan(i/2)*e/(eta^2*G), -m/G/(1+cos(i)), 0, 0, 0];
            Fformx = dFdX;
        else
            Wghcoeffx = [-1/(a*G*4*cos(i/2))*(cos(i)*(n-2*p)-m), e/(eta^2*G*2*cos(i/2))*(cos(i)*(n-2*p)-m), ...
                                    -1/(2*G)*((n-2*p) + (m+n-2*p)/(1+cos(i)))*sin(i/2), 0, 0, 0];
            Fformx = [0,0,Fbysi2d,0,0,0];
        end
        
        Vghx = @(X) CS1*(sbeta*betax(X)*I(2) - cbeta*I2x(X) - cbeta*betax(X)*I(1) - sbeta*I1x(X)) ...
                    + CS2*(-sbeta*betax(X)*I(1) + cbeta*I1x(X) - cbeta*betax(X)*I(2) - sbeta*I2x(X));
        for ctr = 1:6
            Wghx(ctr) = Wghx(ctr) + dCdX(ctr)*Wghcoeff*Fform*Vgh + C*Wghcoeffx(ctr)*Fform*Vgh ...
                                  + C*Wghcoeff*Fformx(ctr)*Vgh + C*Wghcoeff*Fform*Vghx(ctr);
        end
        
        % 1/e*(dWdl - 1/eta*dWdg) term
        Vlgx = @(X) CS1*(-sbeta*betax(X)*(-(n-2*p)*dIdlg1(2) - (n+1)/eta^3*dIdlg2(1)) ...
                            + cbeta*(-(n-2*p)*dI2lgadX(X) - (n+1)*(3*e/eta^5*dIdlg2(1) + 1/eta^3*dI1lgbdX(X))) ...
                         + cbeta*betax(X)*(-(n-2*p)*dIdlg1(1) + (n+1)/eta^3*dIdlg2(2)) ...
                            + sbeta*(-(n-2*p)*dI1lgadX(X) + (n+1)*(3*e/eta^5*dIdlg2(2) + 1/eta^3*dI2lgbdX(X)))) ...
                  + CS2*(-sbeta*betax(X)*((n-2*p)*dIdlg1(1) - (n+1)/eta^3*dIdlg2(2)) ...
                            + cbeta*((n-2*p)*dI1lgadX(X) - (n+1)*(3*e/eta^5*dIdlg2(2) + 1/eta^3*dI2lgbdX(X))) ...
                         + cbeta*betax(X)*(-(n-2*p)*dIdlg1(2) - (n+1)/eta^3*dIdlg2(1)) ...
                            + sbeta*(-(n-2*p)*dI2lgadX(X) - (n+1)*(3*e/eta^5*dIdlg2(1) + 1/eta^3*dI1lgbdX(X))));      

        for ctr = 1:6
            Wlgx(ctr) = Wlgx(ctr) + dCdX(ctr)*F*Vlg + C*dFdX(ctr)*Vlg + C*F*Vlgx(ctr);
        end
    end
    
end
        

end

% function [ J ] = FiniteDiffIxx(X, cdtol, f0,n,m,p,delta,ddeltada,d2deltada2,tol, HessianOn)
% 
% J = zeros(8,6);
% 
% [~, f1] = KeplerEqSolver(X(6), X(2), 1e-11);
% [~, f2] = KeplerEqSolver(X(6) + cdtol, X(2), 1e-11);
% 
% for ctr = 1:6
% 
%     dx1 = X;
%     dx2 = X;
%     dx1(ctr) = X(ctr);
%     dx2(ctr) = X(ctr) + cdtol;
%     
%     dx1(6) = f1;
%     dx2(6) = f1;
%     
%     if ctr == 6
%         dx2(6) = f2;
%     end
% 
%     [I, dIdX, dIdlg1, dIdlg2,I1xx,I2xx,dI1lgdX,dI2lgdX] = qCulpritIntegrals([dx1(1:2),dx1(6)],f0,n,m,p,delta,ddeltada,d2deltada2,tol, HessianOn);
%     dy1 = [dIdX(1,:), dIdlg1(1), dIdlg2(1)]';
%     [I, dIdX, dIdlg1, dIdlg2,I1xx,I2xx,dI1lgdX,dI2lgdX] = qCulpritIntegrals([dx2(1:2),dx2(6)],f0,n,m,p,delta,ddeltada,d2deltada2,tol, HessianOn);
%     dy2 = [dIdX(1,:), dIdlg1(1), dIdlg2(1)]';
%     
%     J(:,ctr) = (dy2 - dy1)/cdtol;
% 
% end
% 
% end
