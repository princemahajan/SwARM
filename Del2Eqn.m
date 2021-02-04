%% Delaunay elemets to Equinoctial Converter
%   Del: [l,g,h,L,G,H];
%   Poin: [a,psi,c(=h),d(=k),p,q];
% Author: Bharat Mahajan (https://github.com/princemahajan)

function [Y] = Del2Eqn(X,mu, toDel)

if toDel == true
   
    a = X(1); Lambda = X(2); p1 = X(3); p2 = X(4); q1 = X(5); q2 = X(6);
    
    % longitude of perihelion
    gh = mod(atan2(q2,q1),2*pi); 
    
    l = Lambda - gh;

    h = mod(atan2(p2,p1),2*pi);
    
    g = mod(gh - h,2*pi);

    L = sqrt(mu*a);
    
    e = sqrt(q1^2 + q2^2);
    
    G = L*sqrt(1 - e^2);
    
    i = 2*atan(sqrt(p1^2 + p2^2));
    
    H = G*cos(i);
    
    Y = [l,g,h,L,G,H]';
    
else
    
    % convert to Eqinoctial elements
    l = X(1); g = X(2); h = X(3); L = X(4); G = X(5); H = X(6);
    
    a = L^2/mu;
    
    Lambda = l + mod(g + h,2*pi);
    
    e = sqrt(1 - G^2/L^2);
    
    q1 = e*cos(g + h);
    
    q2 = e*sin(g + h);
    
    i = acos(H/G);
        
    p1 = tan(i/2)*cos(h);
    p2 = tan(i/2)*sin(h);
    
    Y = [a,Lambda,p1,p2,q1,q2]';
    
end

end
