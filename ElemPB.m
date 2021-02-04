function [ dE, PBJ ] = ElemPB( dWdE,d2WdE2, E, ElemType, mu, JacobianOn)
%ElemPB Evaluates Poisson Brackets of the specified Orbital Elements (functions
%of the canonical Delaunay elements) with respect to a second quantity using
% its partials provided. The second quantity is typically a generating function,
% in that case the output is the variations of the elements due to that
% generating function.
%
% Author: Bharat Mahajan (https://github.com/princemahajan)

if strcmp(ElemType,'classical') == true
    
    % Elements
    a = E(1);
    e = E(2);
    i = E(3);
%     h = E(4);
%     g = E(5);
%     f = E(6);
    
    % partials of the generating function
    dWda = dWdE(1);
    dWde = dWdE(2);
    dWdi = dWdE(3);
    dWdh = dWdE(4);
    dWdg = dWdE(5);
    dWdl = dWdE(6);
    
    % intermediate quantities
    eta = sqrt(1-e^2);
    L = sqrt(mu*a);
    G = L*eta;
    H = G*cos(i);
    
    dadL = 2*L/mu;
    dedL = eta^2/(e*L);
    dedG = -eta^2/(e*G);
    didG = H/(G^2*sin(i));
    didH = -1/(G*sin(i));
    
    % Element Variations
    
    % SMA
    da = -dadL*dWdl;
    
    % Ecc
    de = -dedL*dWdl - dedG*dWdg;
    
    % Inc
    di = -didG*dWdg - didH*dWdh;
    
    % RAAN
    dh = didH*dWdi;
    
    %AOP
    dg = didG*dWdi + dedG*dWde;
    
    % MA
    dl = dWda*dadL + dWde*dedL;
    
    dE = [da, de, di dh, dg, dl]';
    
elseif strcmp(ElemType,'equinoctial') == true

    % Elements
    a = E(1);
    e = E(2);
    i = E(3);
    h = E(4);
    g = E(5);
%     f = E(6);
    
    % partials of the generating function
    dWda = dWdE(1);
    dWde = dWdE(2);
    dWdi = dWdE(3);
    didGWgdidHWh = dWdE(4);
    dWdletadWdg = dWdE(5);
    dWdl = dWdE(6);

    % intermediate quantities
    eta = sqrt(1-e^2);
    L = sqrt(mu*a);
    G = L*eta;
    
    dadL = 2*L/mu;
%     dfdl = (1 + e*cos(f))^2/eta^3;
    
    % Element Variations
    
    % SMA corrections
    DelSMA = -dadL*dWdl;
    
    % Mean Argument of Latitude
    DelMAOL = dadL*dWda - eta/L*dWde*sqrt((1-eta)/(1+eta)) - dWdi*tan(i/2)/G;
    
    % p1 and p2
    Delp1 = sin(h)/(G*(1+cos(i)))*dWdi - cos(h)/(1+cos(i))*didGWgdidHWh;
    Delp2 = -cos(h)/(G*(1+cos(i)))*dWdi - sin(h)/(1+cos(i))*didGWgdidHWh;
    
    % q1,q2
    Delq1 = -sin(g+h)*(-eta/L*dWde-e*dWdi*tan(i/2)/G) - cos(g+h)*eta^2/L*dWdletadWdg;
    Delq2 = cos(g+h)*(-eta/L*dWde-e*dWdi*tan(i/2)/G) - sin(g+h)*eta^2/L*dWdletadWdg;
    
    dE = [DelSMA, DelMAOL, Delp1, Delp2, Delq1, Delq2]';
    
    if JacobianOn == true
        
        % equinoctial and classical element relations
        si = sin(i);
        ci = cos(i);
        ci2 = cos(i/2);
        ti2 = tan(i/2);
        sh = sin(h);
        ch = cos(h);
        cgh = cos(g+h);
        sgh = sin(g+h);
        dedq1 = cgh;
        dedq2 = sgh;
        didp1 = ch*(1+ci);
        didp2 = sh*(1+ci);
        dgdp1 = sh/ti2;
        dgdp2 = -ch/ti2;
        dgdq1 = -sgh/e;
        dgdq2 = cgh/e;
        dhdp1 = -sh/ti2;
        dhdp2 = ch/ti2;
        dldq1 = sgh/e;
        dldq2 = -cgh/e;
        
%         dfde = (2+e*cos(f)).*sin(f)/eta^2;
%         d2fdfdl = 2*(1+e*cos(f))/eta^3*(-e*sin(f));
%         d2fdl2 = dfdl*d2fdfdl;
%         d2fdlde = -3/eta^4*(-e/eta)*(1+e*cos(f))^2 + 2*(1+e*cos(f))/eta^3*(cos(f)-e*sin(f)*dfde);
        
        % Del SMA partials w.r.t. equinoctial elements
        DelSMAa = -1/L*dWdl - dadL*d2WdE2(6,1);
        DelSMAML = -dadL*d2WdE2(6,6);
        DelSMAp1 = -dadL*(d2WdE2(6,3)*didp1 + d2WdE2(6,5)*dgdp1 + d2WdE2(6,4)*dhdp1);
        DelSMAp2 = -dadL*(d2WdE2(6,3)*didp2 + d2WdE2(6,5)*dgdp2 + d2WdE2(6,4)*dhdp2);
        DelSMAq1 = -dadL*((d2WdE2(6,2))*dedq1 + d2WdE2(6,5)*dgdq1 + (d2WdE2(6,6))*dldq1);
        DelSMAq2 = -dadL*((d2WdE2(6,2))*dedq2 + d2WdE2(6,5)*dgdq2 + (d2WdE2(6,6))*dldq2);
    
        % Del MOL partials w.r.t. equinoctial elements
        DelMAOLa = (1/L*dWda+dadL*d2WdE2(1,1)) - eta*sqrt((1-eta)/(1+eta))*(-mu/(2*L^3)*dWde + 1/L*d2WdE2(2,1)) ...
                                - tan(i/2)/eta*(d2WdE2(3,1)/L - dWdi*mu/(2*L^3));
        DelMAOLML = dadL*d2WdE2(1,6) - eta/L*d2WdE2(2,6)*sqrt((1-eta)/(1+eta)) - d2WdE2(3,6)*tan(i/2)/G;
        DelMAOLp1 = dadL*(d2WdE2(1,3)*didp1 + d2WdE2(1,5)*dgdp1 + d2WdE2(1,4)*dhdp1) ...
                    - eta/L*sqrt((1-eta)/(1+eta))*(d2WdE2(2,3)*didp1 + d2WdE2(2,5)*dgdp1 + d2WdE2(2,4)*dhdp1) ...
                    - 1/G*((d2WdE2(3,5)*dgdp1 + d2WdE2(3,4)*dhdp1 + d2WdE2(3,3)*didp1)*tan(i/2) + dWdi*(sec(i/2)^2/2)*didp1);
        DelMAOLp2 = dadL*(d2WdE2(1,3)*didp2 + d2WdE2(1,5)*dgdp2 + d2WdE2(1,4)*dhdp2) ...
                    - eta/L*sqrt((1-eta)/(1+eta))*(d2WdE2(2,3)*didp2 + d2WdE2(2,5)*dgdp2 + d2WdE2(2,4)*dhdp2) ...
                    - 1/G*((d2WdE2(3,5)*dgdp2 + d2WdE2(3,4)*dhdp2 + d2WdE2(3,3)*didp2)*tan(i/2) + dWdi*(sec(i/2)^2/2)*didp2);
        DelMAOLq1 = dadL*(d2WdE2(1,2)*dedq1 + d2WdE2(1,5)*dgdq1 + d2WdE2(1,6)*dldq1)...
                    - 1/L*((eta*sqrt((1-eta)/(1+eta))*d2WdE2(2,2) + ((2*eta-1)/eta-eta/(1+eta))*dWde)*dedq1 ...
                            + eta*sqrt((1-eta)/(1+eta))*(d2WdE2(2,5)*dgdq1 + d2WdE2(2,6)*dldq1))  ...
                    - 1/L*((-1/eta^2)*(-e/eta)*dedq1*tan(i/2)*dWdi + 1/eta*tan(i/2)*d2WdE2(3,2)*dedq1 ...
                            + 1/eta*tan(i/2)*(d2WdE2(3,5)*dgdq1 + d2WdE2(3,6)*dldq1));
        DelMAOLq2 = dadL*(d2WdE2(1,2)*dedq2 + d2WdE2(1,5)*dgdq2 + d2WdE2(1,6)*dldq2)...
                    - 1/L*((eta*sqrt((1-eta)/(1+eta))*d2WdE2(2,2) + ((2*eta-1)/eta-eta/(1+eta))*dWde)*dedq2 ...
                            + eta*sqrt((1-eta)/(1+eta))*(d2WdE2(2,5)*dgdq2 + d2WdE2(2,6)*dldq2))  ...
                    - 1/L*((-1/eta^2)*(-e/eta)*dedq2*tan(i/2)*dWdi + 1/eta*tan(i/2)*d2WdE2(3,2)*dedq2 ...
                            + 1/eta*tan(i/2)*(d2WdE2(3,5)*dgdq2 + d2WdE2(3,6)*dldq2));
    
        % Del p1 and p2 partials w.r.t. equinoctial elements
        Delp1a = sh/(1+ci)/eta*(-mu/(2*L^3)*dWdi + 1/L*d2WdE2(3,1)) - ch/(1+ci)*d2WdE2(7,1);
        Delp2a = -ch/(1+ci)/eta*(-mu/(2*L^3)*dWdi + 1/L*d2WdE2(3,1)) - sh/(1+ci)*d2WdE2(7,1);
        
        Delp1ML = sh/(G*(1+ci))*d2WdE2(3,6) - ch/(1+ci)*d2WdE2(7,6);
        Delp2ML = -ch/(G*(1+ci))*d2WdE2(3,6) - sh/(1+ci)*d2WdE2(7,6);
        
        Delp1p1 = 1/G*(sh*(si/(1+ci)^2*dWdi + 1/(1+ci)*d2WdE2(3,3))*didp1 + 1/(1+ci)*(sh*d2WdE2(3,5)*dgdp1 + (ch*dWdi + sh*d2WdE2(3,4))*dhdp1)) ...
                 -ch*(si/(1+ci)^2*didGWgdidHWh + 1/(1+ci)*d2WdE2(7,3))*didp1 - 1/(1+ci)*(ch*d2WdE2(7,5)*dgdp1+(-sh*didGWgdidHWh+ch*d2WdE2(7,4))*dhdp1);
        Delp1p2 = 1/G*(sh*(si/(1+ci)^2*dWdi + 1/(1+ci)*d2WdE2(3,3))*didp2 + 1/(1+ci)*(sh*d2WdE2(3,5)*dgdp2 + (ch*dWdi + sh*d2WdE2(3,4))*dhdp2)) ...
                 -ch*(si/(1+ci)^2*didGWgdidHWh + 1/(1+ci)*d2WdE2(7,3))*didp2-1/(1+ci)*(ch*d2WdE2(7,5)*dgdp2+(-sh*didGWgdidHWh+ch*d2WdE2(7,4))*dhdp2);
        
        Delp2p1 = 1/G*(-ch*(si/(1+ci)^2*dWdi + 1/(1+ci)*d2WdE2(3,3))*didp1 + 1/(1+ci)*(-ch*d2WdE2(3,5)*dgdp1 + (sh*dWdi - ch*d2WdE2(3,4))*dhdp1)) ...
                 -sh*(si/(1+ci)^2*didGWgdidHWh + 1/(1+ci)*d2WdE2(7,3))*didp1 - 1/(1+ci)*(sh*d2WdE2(7,5)*dgdp1+(ch*didGWgdidHWh+sh*d2WdE2(7,4))*dhdp1);
        Delp2p2 = 1/G*(-ch*(si/(1+ci)^2*dWdi + 1/(1+ci)*d2WdE2(3,3))*didp2 + 1/(1+ci)*(-ch*d2WdE2(3,5)*dgdp2 + (sh*dWdi - ch*d2WdE2(3,4))*dhdp2)) ...
                 -sh*(si/(1+ci)^2*didGWgdidHWh + 1/(1+ci)*d2WdE2(7,3))*didp2 - 1/(1+ci)*(sh*d2WdE2(7,5)*dgdp2+(ch*didGWgdidHWh+sh*d2WdE2(7,4))*dhdp2);
        
        Delp1q1 = 1/L*sh/(1+ci)*((-1/eta^2*(-e/eta)*dWdi + 1/eta*d2WdE2(3,2))*dedq1 + 1/eta*(d2WdE2(3,5)*dgdq1 + d2WdE2(3,6)*dldq1)) ...
                    - ch/(1+ci)*(d2WdE2(7,2)*dedq1 + d2WdE2(7,5)*dgdq1 + d2WdE2(7,6)*dldq1);
        Delp1q2 = 1/L*sh/(1+ci)*((-1/eta^2*(-e/eta)*dWdi + 1/eta*d2WdE2(3,2))*dedq2 + 1/eta*(d2WdE2(3,5)*dgdq2 + d2WdE2(3,6)*dldq2)) ...
                    - ch/(1+ci)*(d2WdE2(7,2)*dedq2 + d2WdE2(7,5)*dgdq2 + d2WdE2(7,6)*dldq2);
        
        Delp2q1 = -1/L*ch/(1+ci)*((-1/eta^2*(-e/eta)*dWdi + 1/eta*d2WdE2(3,2))*dedq1 + 1/eta*(d2WdE2(3,5)*dgdq1 + d2WdE2(3,6)*dldq1)) ...
                    - sh/(1+ci)*(d2WdE2(7,2)*dedq1 + d2WdE2(7,5)*dgdq1 + d2WdE2(7,6)*dldq1);
        Delp2q2 = -1/L*ch/(1+ci)*((-1/eta^2*(-e/eta)*dWdi + 1/eta*d2WdE2(3,2))*dedq2 + 1/eta*(d2WdE2(3,5)*dgdq2 + d2WdE2(3,6)*dldq2)) ...
                    - sh/(1+ci)*(d2WdE2(7,2)*dedq2 + d2WdE2(7,5)*dgdq2 + d2WdE2(7,6)*dldq2);
        

        % Del q1,q2 Partials w.r.t equinoctial elements
        Delq1a = -sgh*((-mu/(2*L^3)*(-eta*dWde-e*dWdi*ti2/eta) + 1/L*(-eta*d2WdE2(2,1)-e*d2WdE2(3,1)*ti2/eta))) ...
                        -cgh*eta^2*(-mu/(2*L^3)*dWdletadWdg + 1/L*d2WdE2(8,1));
        Delq2a = cgh*((-mu/(2*L^3)*(-eta*dWde-e*dWdi*ti2/eta) + 1/L*(-eta*d2WdE2(2,1)-e*d2WdE2(3,1)*ti2/eta))) ...
                        -sgh*eta^2*(-mu/(2*L^3)*dWdletadWdg + 1/L*d2WdE2(8,1));
                    
        Delq1ML = -sgh*(-eta/L*d2WdE2(2,6)-e*d2WdE2(3,6)*tan(i/2)/G) - cgh*eta^2/L*d2WdE2(8,6);
        Delq2ML = cgh*(-eta/L*d2WdE2(2,6)-e*d2WdE2(3,6)*tan(i/2)/G) - sgh*eta^2/L*d2WdE2(8,6);
        
        Delq1p1 = -sgh*(-eta/L*d2WdE2(2,3)-e*d2WdE2(3,3)*ti2/G-e*dWdi/(2*ci2^2*G))*didp1 - cgh*(dgdp1+dhdp1)*(-eta/L*dWde-e*dWdi*ti2/G) ...
                            -sgh*(-eta/L*(d2WdE2(2,5)*dgdp1+d2WdE2(2,4)*dhdp1)-e*ti2/G*(d2WdE2(3,5)*dgdp1+d2WdE2(3,4)*dhdp1)) ...
                    -eta^2/L*(cgh*d2WdE2(8,3)*didp1  -sgh*(dgdp1+dhdp1)*dWdletadWdg + cgh*(d2WdE2(8,5)*dgdp1+d2WdE2(8,4)*dhdp1));
        Delq2p1 = cgh*(-eta/L*d2WdE2(2,3)-e*d2WdE2(3,3)*ti2/G-e*dWdi/(2*ci2^2*G))*didp1 - sgh*(dgdp1+dhdp1)*(-eta/L*dWde-e*dWdi*ti2/G) ...
                            +cgh*(-eta/L*(d2WdE2(2,5)*dgdp1+d2WdE2(2,4)*dhdp1)-e*ti2/G*(d2WdE2(3,5)*dgdp1+d2WdE2(3,4)*dhdp1)) ...
                    -eta^2/L*(sgh*d2WdE2(8,3)*didp1 + cgh*(dgdp1+dhdp1)*dWdletadWdg + sgh*(d2WdE2(8,5)*dgdp1+d2WdE2(8,4)*dhdp1));
        
        Delq1p2 = -sgh*(-eta/L*d2WdE2(2,3)-e*d2WdE2(3,3)*ti2/G-e*dWdi/(2*ci2^2*G))*didp2 - cgh*(dgdp2+dhdp2)*(-eta/L*dWde-e*dWdi*ti2/G) ...
                            -sgh*(-eta/L*(d2WdE2(2,5)*dgdp2+d2WdE2(2,4)*dhdp2)-e*ti2/G*(d2WdE2(3,5)*dgdp2+d2WdE2(3,4)*dhdp2)) ...
                    -eta^2/L*(cgh*d2WdE2(8,3)*didp2 -sgh*(dgdp2+dhdp2)*dWdletadWdg + cgh*(d2WdE2(8,5)*dgdp2+d2WdE2(8,4)*dhdp2));
        Delq2p2 = cgh*(-eta/L*d2WdE2(2,3)-e*d2WdE2(3,3)*ti2/G-e*dWdi/(2*ci2^2*G))*didp2 - sgh*(dgdp2+dhdp2)*(-eta/L*dWde-e*dWdi*ti2/G) ...
                            +cgh*(-eta/L*(d2WdE2(2,5)*dgdp2+d2WdE2(2,4)*dhdp2)-e*ti2/G*(d2WdE2(3,5)*dgdp2+d2WdE2(3,4)*dhdp2)) ...
                    -eta^2/L*(sgh*d2WdE2(8,3)*didp2 +cgh*(dgdp2+dhdp2)*dWdletadWdg + sgh*(d2WdE2(8,5)*dgdp2+d2WdE2(8,4)*dhdp2));
        
        Delq1q1 = -1/L*(sgh*(e/eta*dWde-eta*d2WdE2(2,2)-(1/eta^3*dWdi+e/eta*d2WdE2(3,2))*ti2)*dedq1 ...
                        +(cgh*dgdq1*(-eta*dWde-e*dWdi*ti2/eta)+sgh*((-eta*d2WdE2(2,5)-e*d2WdE2(3,5)*ti2/eta)*dgdq1...
                        +(-eta*d2WdE2(2,6)-e*d2WdE2(3,6)*ti2/eta)*dldq1))) ...
                    - 1/L*((cgh*2*eta*(-e/eta)*dWdletadWdg + cgh*eta^2*d2WdE2(8,2))*dedq1 + eta^2*(-sgh*dWdletadWdg + cgh*d2WdE2(8,5))*dgdq1 ...
                           + eta^2*cgh*d2WdE2(8,6)*dldq1);
        Delq2q1 = -1/L*(-cgh*(e/eta*dWde-eta*d2WdE2(2,2)-(1/eta^3*dWdi+e/eta*d2WdE2(3,2))*ti2)*dedq1 ...
                        +(sgh*dgdq1*(-eta*dWde-e*dWdi*ti2/eta)-cgh*((-eta*d2WdE2(2,5)-e*d2WdE2(3,5)*ti2/eta)*dgdq1...
                        +(-eta*d2WdE2(2,6)-e*d2WdE2(3,6)*ti2/eta)*dldq1))) ...
                    - 1/L*((sgh*2*eta*(-e/eta)*dWdletadWdg + sgh*eta^2*d2WdE2(8,2))*dedq1 + eta^2*(cgh*dWdletadWdg + sgh*d2WdE2(8,5))*dgdq1 ...
                           + eta^2*sgh*d2WdE2(8,6)*dldq1);
        
         Delq1q2 = -1/L*(sgh*(e/eta*dWde-eta*d2WdE2(2,2)-(1/eta^3*dWdi+e/eta*d2WdE2(3,2))*ti2)*dedq2 ...
                        +(cgh*dgdq2*(-eta*dWde-e*dWdi*ti2/eta)+sgh*((-eta*d2WdE2(2,5)-e*d2WdE2(3,5)*ti2/eta)*dgdq2...
                        +(-eta*d2WdE2(2,6)-e*d2WdE2(3,6)*ti2/eta)*dldq2))) ...
                    - 1/L*((cgh*2*eta*(-e/eta)*dWdletadWdg + cgh*eta^2*d2WdE2(8,2))*dedq2 + eta^2*(-sgh*dWdletadWdg + cgh*d2WdE2(8,5))*dgdq2 ...
                           + eta^2*cgh*d2WdE2(8,6)*dldq2);
        Delq2q2 = -1/L*(-cgh*(e/eta*dWde-eta*d2WdE2(2,2)-(1/eta^3*dWdi+e/eta*d2WdE2(3,2))*ti2)*dedq2 ...
                        +(sgh*dgdq2*(-eta*dWde-e*dWdi*ti2/eta)-cgh*((-eta*d2WdE2(2,5)-e*d2WdE2(3,5)*ti2/eta)*dgdq2...
                        +(-eta*d2WdE2(2,6)-e*d2WdE2(3,6)*ti2/eta)*dldq2))) ...
                    - 1/L*((sgh*2*eta*(-e/eta)*dWdletadWdg + sgh*eta^2*d2WdE2(8,2))*dedq2 + eta^2*(cgh*dWdletadWdg + sgh*d2WdE2(8,5))*dgdq2 ...
                           + eta^2*sgh*d2WdE2(8,6)*dldq2);
        
        
        % Short-period Jacobian
        PBJ = [DelSMAa, DelSMAML, DelSMAp1, DelSMAp2, DelSMAq1, DelSMAq2;
               DelMAOLa, DelMAOLML, DelMAOLp1, DelMAOLp2, DelMAOLq1, DelMAOLq2;
               Delp1a, Delp1ML, Delp1p1, Delp1p2, Delp1q1, Delp1q2;
               Delp2a, Delp2ML, Delp2p1, Delp2p2, Delp2q1, Delp2q2;
               Delq1a, Delq1ML, Delq1p1, Delq1p2, Delq1q1, Delq1q2;
               Delq2a, Delq2ML, Delq2p1, Delq2p2, Delq2q1, Delq2q2];
        
    else
        PBJ = 0;
    end
    
end

