%% ECI2LVLH - ECI to Curvilinear frame conversion and vice-versa
%   Inputs:
%       X1: Reference Frame Pos,Vel in ECI
%       X2: Object Pos,Vel in ECI if reverse not true
%           Object relative Pos,Vel in Curvilinaer frame if reverse true
%       reverse: true for converting from LVLH to ECI
%   Units: rad, km, sec 
%   Tolerance: 1e-12
% Author: Bharat Mahajan (https://github.com/princemahajan)

function [ Y ] = ECI2Curv( X1, X2, reverse )

if size(X2,2) ~= 1
    disp('ECI2Curv: Second parameter must be a vector/column matrix');
    return;
end

% form DCM
rcap = X1(1:3,1)'/norm(X1(1:3));
hvec = cross(X1(1:3,1)',X1(4:6,1)');
hcap = hvec/norm(hvec);
tcap = cross(hcap,rcap);
DCM = [rcap; tcap; hcap];

% Angular velocity In Hill frame of the chief frame (r^2*w = h)
w = [0;0;norm(hvec)/norm(X1(1:3))^2];

if reverse == true
   
    % convert from Curvilinear to Hill's frame
    XLVLH = curv2Hill(X2, X1);
    
    % ECI pos
    Y(1:3,1) = DCM'*XLVLH(1:3,1) + X1(1:3,1);
    
    % Inertial relative velocity expressed in LVLH
    Y(4:6,1) = XLVLH(4:6,1) + cross(w,XLVLH(1:3,1)); 
    
    % ECI velocity in ECI frame
    Y(4:6,1) = DCM'*Y(4:6,1) + X1(4:6,1);
    
else
    
    % convert from ECI to LVLH
    
    % inertial difference vector
    dX = X2 - X1;
    
    % LVLH pos
    Y(1:3,1) = DCM*dX(1:3,1);
    
    % LVLH vel
    Y(4:6,1) = DCM*dX(4:6,1) - cross(w,Y(1:3,1));
   
    % to curvilinear frame
    Y = Hill2curv(Y, X1);
end



end

