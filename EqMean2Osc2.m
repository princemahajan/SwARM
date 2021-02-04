function [ Y,Ylp,D,Dlp] = EqMean2Osc2( X, GMST, mu, Re, we, degree, order, Clm, Slm, tol, quadtol, InverseOn, UseMex, JacobianOn, SEMI_ANALYTIC, QuadTesseralsOn )
%Mean2Osc Mean To Osculating Transformation that includes J2^2 short-period and
%long period, other zonals and tesserals short-period up to 2nd order and 
%long period up to 1st order. 
%
% Uses Equinoctial elements:
% [a,l+g+h,tan(i/2)cos(h),tan(i/2)sin(h),ecos(g+h),esin(g+h)]
%
% Author: Bharat Mahajan (https://github.com/princemahajan)

% IF Two-body propagation
if degree < 2
   Y = X;
   Ylp = X;
   D = eye(6);
   Dlp = eye(6);
   return;
end

% Tesseral File
% EqSPMDTesseralsFile = @EqSPMDTesseralsMRe;
EqSPMDTesseralsFile = @EqSPMDTesserals;

% Coefficients
Jcoeff = -Clm(2:end,:);
J2 = Jcoeff(2);

% Choose whether use MEX files for J2 LP and SP Code
if UseMex == true
    EqJ2LPTransform = @EqLPJ2Mex;
    EqJ2SPTransform = @EqSPJ2Mex;
else
    EqJ2LPTransform = @EqLPJ2;
    EqJ2SPTransform = @EqSPJ2;
end

% Decide which transform-LP or SP to apply first?
if InverseOn == true
    J2Trans1 = EqJ2SPTransform;
    if SEMI_ANALYTIC == true
        J2Trans2 = @dummy1;
    else
        J2Trans2 = EqJ2LPTransform;
    end
    JnTrans1 = @EqSPZonal;
    
    if SEMI_ANALYTIC == true
        JnTrans2 = @dummy2;
    else
        JnTrans2 = @EqLPZonal;
    end

    % Equinoctial m-Daily and Tesseral SP Variations
    TessX = X;
    [DelTessSP, DspT] = EqSPMDTesseralsFile(TessX, GMST, we, mu, Re, Clm, Slm, degree, order, tol, quadtol,JacobianOn,QuadTesseralsOn);
    % Tesseral D matrix
%     cdtol = 1e-10;
    if JacobianOn == true
%         DspT = FiniteDiff(TessX, EqSPMDTesseralsFile, cdtol, GMST, we, mu, Re, Clm, Slm, degree, order, tol, quadtol );
%         DspT(1,5:6) = DspT1(1,5:6);
%         DspT(5,:) = DspT1(5,:);
%         DspT(6,:) = DspT1(6,:);
    else
        DspT = 0;
    end
    
    zX = X - DelTessSP;
    
    % Equinoctial Zonal 1st Variations
    [Del1J2, D1J2] = J2Trans1(zX, mu, Re, J2, InverseOn, JacobianOn, tol);
    [Del1Jn, D1Jn] = JnTrans1(zX, degree ,mu, Re, Jcoeff, InverseOn, JacobianOn, tol);
    
    Del1 = Del1J2 + Del1Jn;
    D1 = D1J2 + D1Jn;
    
    % Secular + LP states (Zonal)
    Ylp = zX + Del1;
    Dlp = eye(6) + D1;
    
    % Equinoctial Zonal 2nd Variations
    [Del2J2, D2J2] = J2Trans2(Ylp, mu,Re,J2, InverseOn, UseMex,tol);
    [Del2Jn, D2Jn] = JnTrans2(Ylp, degree,mu,Re,Jcoeff, InverseOn, UseMex, tol);
    
    Del2 = Del2J2 + Del2Jn;
    D2 = D2J2 + D2Jn;
    
    % Add all the variations
    Y = zX + Del1 + Del2;
    
    % Jacobian
    D = (eye(6) + D2)*(eye(6) + D1)*(eye(6) + DspT)^-1;
    
    
else
    
    if SEMI_ANALYTIC == true
        J2Trans1 = @dummy1;
    else
        J2Trans1 = EqJ2LPTransform;
    end
    J2Trans2 = EqJ2SPTransform;

    if SEMI_ANALYTIC == true
        JnTrans1 = @dummy2;
    else
        JnTrans1 = @EqLPZonal;
    end
    JnTrans2 = @EqSPZonal;
    
    % Equinoctial Zonal 1st Variations
    [Del1J2, D1J2] = J2Trans1(X, mu, Re, J2, InverseOn, JacobianOn, tol);
    [Del1Jn, D1Jn] = JnTrans1(X, degree ,mu, Re, Jcoeff, InverseOn, JacobianOn, tol);
    
    Del1 = Del1J2 + Del1Jn;
    D1 = D1J2 + D1Jn;
    
    % Secular + LP states (Zonal)
    Ylp = X + Del1;
    Dlp = eye(6) + D1;
    
    % Equinoctial Zonal 2nd Variations
    [Del2J2, D2J2] = J2Trans2(Ylp, mu,Re,J2, InverseOn, UseMex,tol);
    [Del2Jn, D2Jn] = JnTrans2(Ylp, degree,mu,Re,Jcoeff, InverseOn, UseMex, tol);
    
    Del2 = Del2J2 + Del2Jn;
    D2 = D2J2 + D2Jn;
    
    % Equinoctial m-Daily and Tesseral SP Variations
    TessX = X + Del1 + Del2;
    [DelTessSP, DspT] = EqSPMDTesseralsFile(TessX, GMST, we, mu, Re, Clm, Slm, degree, order, tol, quadtol, JacobianOn,QuadTesseralsOn);
    % Tesseral D matrix
%     cdtol = 1e-10;
    if JacobianOn == true
%         DspT = FiniteDiff(TessX, EqSPMDTesseralsFile, cdtol, GMST, we, mu, Re, Clm, Slm, degree, order, tol, quadtol );
%         DspT(1,5:6) = DspT1(1,5:6);
%         DspT(5,:) = DspT1(5,:);
%         DspT(6,:) = DspT1(6,:);
    else
        DspT = 0;
    end
    
    % Add all the variations
    Y = X + Del1 + Del2 + DelTessSP;
    
    % Jacobian
    D = (eye(6) + DspT)*(eye(6) + D2)*(eye(6) + D1);
    
end


end

% function [ J ] = FiniteDiff(X, VecFunc, cdtol, GMST, we, mu, Re, Clm, Slm, degree, order, tol, quadtol )
% 
% J = zeros(6);
% 
% for ctr = 1:6
% 
%     dx1 = X;
%     dx2 = X;
%     dx1(ctr) = X(ctr);
%     dx2(ctr) = X(ctr) + cdtol;
% 
%     dy1 = VecFunc(dx1, GMST, we, mu, Re, Clm, Slm, degree, order, tol, quadtol, false );
%     dy2 = VecFunc(dx2, GMST, we, mu, Re, Clm, Slm, degree, order, tol, quadtol, false );
%     
%     J(:,ctr) = (dy2 - dy1)/cdtol;
% 
% end
% 
% end
% 
% function [Del1J2, D1J2] = dummy1(Ylp, mu,Re,J2, InverseOn, UseMex,tol)
% 
% Del1J2 = 0;
% D1J2 = 0;
% 
% end
% 
% function [Del2Jn, D2Jn] = dummy2(X, degree ,mu, Re, Jcoeff, InverseOn, JacobianOn, tol);
% 
% Del2Jn = 0;
% D2Jn = 0;
% 
% end
% 
