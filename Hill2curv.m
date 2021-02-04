function x_c = Hill2curv(x_h, Xc)
% x_c = Hill2curv(x_h, Xc)
%
% Converts deputy state in Hill frame to curvilinear frame
%
% x_h = [r_h, v_h] -> deputy state in Hill frame
% x_c = [r_c, v_c] -> deputy state in curvilinear frame
% Xc = [R, V] -> chief state in inertial frame
%
% Author: Bharat Mahajan (https://github.com/princemahajan)
%

Rc = norm(Xc(1:3));
Rd_vec=[Rc+x_h(1);x_h(2);x_h(3)];

Rd=norm(Rd_vec);
Rcd = (Xc(1:3)'*Xc(4:6))/Rc;
Rd_dot_vec=[Rcd+x_h(4);x_h(5);x_h(6)];
Rdd=Rd_vec'*Rd_dot_vec/Rd;

x_c = zeros(6,1);

% Position
phi = atan2(x_h(2),(Rc+x_h(1)));
%psi = asin(x_h(3)/Rd);
psi=atan2(x_h(3),sqrt((Rc+x_h(1))^2+x_h(2)^2));

x_c(1) = Rd - Rc;
x_c(2) = Rc*phi;
x_c(3) = Rd*psi;

% Velocity
phid = (x_h(5)*(Rc+x_h(1))-x_h(2)*(Rcd+x_h(4)))/((Rc+x_h(1))^2+x_h(2)^2);
psid = (x_h(6)*sqrt((Rc+x_h(1))^2+x_h(2)^2)-x_h(3)*((Rc+x_h(1))^2+x_h(2)^2)^(-.5)*((Rc+x_h(1))*(Rcd+x_h(4))+x_h(2)*x_h(5)))/Rd^2;

x_c(4) = Rdd-Rcd;
x_c(5) = Rcd*phi + Rc*phid;
x_c(6) = Rdd*psi + psid*Rd;