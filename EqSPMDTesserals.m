%% EqSPMDTesserals
%  Desc: Computes 2nd order Short-period m-daily and Tesseral corrections for 
%        equinoctial elements. See Exact Normalization of Tesseral
%        Harmonics by Bharat Mahajan, Srinivas R. Vadali, Kyle T. Alfriend.
% Author: Bharat Mahajan (https://github.com/princemahajan)

function [DelXm, JdX] = EqSPMDTesserals(Xm,GMST,we,mu,Re,Clm,Slm,MaxDegree,MaxOrder,tol,quadtol, JacobianOn,QuadTesseralsOn) 

% Equinoctial mean elements
a = Xm(1); 
% MAOL = Xm(2); 
% p1 = Xm(3); 
% p2 = Xm(4); 
% q1 = Xm(5); 
% q2 = Xm(6);
 
% current Greenwich mean sidereal time
theta = GMST;

% Specified order must be smaller than degree
MaxOrder = min([MaxDegree, MaxOrder]);

% compute Delaunay elements
DXm = Del2Eqn(Xm,mu, true);
l = DXm(1);
g = mod(DXm(2),2*pi);
h = mod(DXm(3),2*pi);
L = DXm(4);
G = DXm(5);
H = DXm(6);

eta = G/L; 
e = sqrt(1 - eta^2);
i = acos(H/G);

[~, f] = KeplerEqSolver(l, e, tol);

f0 = 0;

DelXm = zeros(6,1);
JdX = zeros(6);

% corrections for each m-Daily, sectorial and tesseral
for n = 2:1:MaxDegree
    for m = 1:1:min(MaxOrder,n)
        
        % C and S coefficients
        C = Clm(n + 1,m + 1);
        S = Slm(n + 1,m + 1);
        C20 = Clm(3,1);

        % short-period and m-daily variations for each tesseral
        coe = [a,e,i,h,g,f]';
        hr = h - theta;
        coer = [a,e,i,hr,g,f]';
%         coerl = [a,e,i,hr,g,l]';
        
        % GF Partials in the ECF frame
        [~, Wx, didGWgdidHWh, dWdletadWdg,Wxx,Wghx,Wlgx] = TesseralQGFr(coer,f0,n,m,C20,C,S,mu,we,Re,quadtol,JacobianOn,QuadTesseralsOn);
        
        % Compute Poisson brackets for equinoctial elements
        dWdE = [Wx(1:3), didGWgdidHWh, dWdletadWdg, Wx(6)];

        % Generating Function Jacobian
        d2WdE2 = [Wxx; Wghx; Wlgx];
%         if JacobianOn == true
%             cdtol = 1e-10;
%             d2WdE2 = FiniteDiffWxx(coerl, cdtol, f0,n,m,C20,C,S,mu,we,Re,quadtol,false );
%         end
        
        [dEqn, PBJ] = ElemPB( dWdE, d2WdE2, coe, 'equinoctial', mu, JacobianOn);
        
        Delnm = C20^2/2*dEqn;
        PBJnm = C20^2/2*PBJ;
        
        DelXm = DelXm + Delnm;
        JdX = JdX + PBJnm;
    end
end

end

% function [ J ] = FiniteDiffWxx(X, cdtol, f0,n,m,C20,C,S,mu,we,Re,quadtol,false )
% 
% J = zeros(8,6);
% 
% [~, f1] = KeplerEqSolver(X(6), X(2), 1e-11);
% [~, f2] = KeplerEqSolver(X(6)+cdtol, X(2), 1e-11);
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
% %     dy1 = VecFunc(dx1, GMST, we, mu, Re, Clm, Slm, degree, order, tol, quadtol, false );
% %     dy2 = VecFunc(dx2, GMST, we, mu, Re, Clm, Slm, degree, order, tol, quadtol, false );
%     
%     [~, Wx, didGWgdidHWh, dWdletadWdg,Wxx,Wghx,Wlgx] = TesseralQGFr(dx1,f0,n,m,C20,C,S,mu,we,Re,quadtol,false);
%     dy1 = [Wx, didGWgdidHWh, dWdletadWdg]';
%     [~, Wx, didGWgdidHWh, dWdletadWdg,Wxx,Wghx,Wlgx] = TesseralQGFr(dx2,f0,n,m,C20,C,S,mu,we,Re,quadtol,false);
%     dy2 = [Wx, didGWgdidHWh, dWdletadWdg]';
%     
%     
%     J(:,ctr) = (dy2 - dy1)/cdtol;
% 
% end
% 
% end
