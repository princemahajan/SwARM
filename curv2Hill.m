function x_h = curv2Hill(xc1, x_eci)
% Dec 25, 2012
% x_h = CURV2HILL(x_c,Rc)
%
% Determines the relative position and velocity of the deputy in the
% LVLH frame using the curvilinear coordinates
%
% Inputs: x_c --> relative position and velocity in the curv frame
%                   [x xdot y ydot z zdot]
%         Rc  --> position of the chief 
% Outputs: x_h --> relative position and velocity in the LVLH frame
%                   [x y z xdot ydot zdot]
% Author: Bharat Mahajan (https://github.com/princemahajan)

xc(1,1) = xc1(1);
xc(2,1) = xc1(4);
xc(3,1) = xc1(2);
xc(4,1) = xc1(5);
xc(5,1) = xc1(3);
xc(6,1) = xc1(6);

R  = x_eci(1:3); V = x_eci(4:6);

x_h = xc;

Rc  = sqrt(R'*R);
Rd  = Rc + xc(1);
phi = xc(3)/Rc;
psi = xc(5)/Rc;

Rc_dot = R'*V/Rc;
Rd_dot = Rc_dot + xc(2);
phidot = (Rc*xc(4) - xc(3)*Rc_dot)/Rc^2;
psidot = (Rc*xc(6) - xc(5)*Rc_dot )/Rc^2;

% Position conversion
x_h(1:3) = [Rd*cos(psi)*cos(phi) - Rc; Rd*cos(psi)*sin(phi);Rd*sin(psi)];

% Velocity conversion
x_h(4) = Rd_dot*cos(phi)*cos(psi) - Rd*sin(psi)*cos(phi)*psidot - Rd*sin(phi)*cos(psi)*phidot - Rc_dot;
x_h(5) = Rd_dot*sin(phi)*cos(psi) - Rd*sin(psi)*sin(phi)*psidot + Rd*cos(phi)*cos(psi)*phidot;
x_h(6) = Rd_dot*sin(psi) + Rd*cos(psi)*psidot;