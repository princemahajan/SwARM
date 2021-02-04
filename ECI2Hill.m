%% ECI2LVLH - ECI to LVLH frame conversion and vice-versa
%   Inputs:
%       X1: Reference Frame Pos,Vel in ECI
%       X2: Object Pos,Vel in ECI if reverse not true
%           Object relative Pos,Vel in LVLH frame if reverse true
%       Omega: Reference frame angular velocity expressed in LVLH frame
%       reverse: true for converting from LVLH to ECI
%   Units: rad, km, sec 
%   Tolerance: 1e-12
% Author: Bharat Mahajan (https://github.com/princemahajan)

function [ Y ] = ECI2Hill( X1, X2, reverse )

% form DCM
rcap = X1(1:3,1)'/norm(X1(1:3));
hvec = cross(X1(1:3,1)',X1(4:6,1)');
hcap = hvec/norm(hvec);
tcap = cross(hcap,rcap);
DCM = [rcap; tcap; hcap];

Omega = [0;0;norm(hvec)/norm(X1(1:3))^2];

if reverse == true
   
    % convert from LVLH to ECI
    
    % ECI pos
    Y(1:3,1) = DCM'*X2(1:3,1) + X1(1:3,1);
    
    % Inertial relative velocity expressed in LVLH
    Y(4:6,1) = X2(4:6,1) + cross(Omega,X2(1:3,1)); 
    
    % ECI velocity in ECI frame
    Y(4:6,1) = DCM'*Y(4:6,1) + X1(4:6,1);
    
else
    
    % convert from ECI to LVLH
    
    % inertial difference vector
    dX = X2 - X1;
    
    % LVLH pos
    Y(1:3,1) = DCM*dX(1:3,1);
    
    % LVLH vel
    Y(4:6,1) = DCM*dX(4:6,1) - cross(Omega,Y(1:3,1));
    
end



end

