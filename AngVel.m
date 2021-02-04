%% AngVel
% AngVel: Computes osculating angular velocity vector in the Hill's frame of a 
% satellite in the presence of zonal and tesseral harmonic perturbations, 
% given satellite's orbital elements.
%
% Author: Bharat Mahajan (https://github.com/princemahajan)

function w = AngVel(OE,theta,degree,order, mu, Re, Cnm,Snm) 

order = min([degree, order]);

a = OE(1);
e = OE(2);
i = OE(3);
h = OE(4);
g = OE(5);
f = OE(6);

% Angular Velocity in x direction with all the zonals and tesserals included
eta = sqrt(1 - e^2);
sv = sin(f + g);
cv = cos(f + g);
si = sin(i);
ci = cos(i);
b = a*sqrt(1 - e^2);
mm = sqrt(mu/a^3);
r = a*eta^2/(1+e*cos(f));

% radial angular velocity: wr = RAANdot*sin(i)/sin(f+g)
wr = 0;

for n = 2:degree
    for m = 0:min(order,n)
     
        % Gravity coefficients
        C = Cnm(n+1,m+1);
        if m == 0
            S = 0;
        else
            S = Snm(n+1,m+1);
        end
        
        CSterm1 = C*cos(m*(h - theta)) + S*sin(m*(h - theta));
        CSterm2 = C*sin(m*(h - theta)) - S*cos(m*(h - theta));
        
        k = fix((n-m)/2);

        wrts = 0;
        
        for t = 0:k
            for s = 0:m
                
                % See Kaula book chapter-1
                T = (-1)^t*factorial(2*n - 2*t)/(2^n*factorial(t)*factorial(n-t)*factorial(n-m-2*t));

                % choose real coefficients
                if mod(s,2) == 0
                    CSterm = (1i^s)*CSterm1;
                else
                    CSterm = (1i^(s+1))*CSterm2;
                end
                
                if (n-m-2*t) == 0 && s == 0
                    wrts = wrts + 0;
                elseif (n-m-2*t) == 0 && s ~= 0
                    wrts = wrts + T*nchoosek(m,s)*CSterm*sv^(n-m-2*t+s-1)*cv^(m-s)*(si^(n-m-2*t+1)*ci^(s-1)*(-s));
                elseif (n-m-2*t) ~= 0 && s == 0
                    wrts = wrts + T*nchoosek(m,s)*CSterm*sv^(n-m-2*t+s-1)*cv^(m-s)*(si^(n-m-2*t-1)*ci^(s+1)*(n-m-2*t));
                else
                    wrts = wrts + T*nchoosek(m,s)*CSterm*sv^(n-m-2*t+s-1)*cv^(m-s)*(si^(n-m*2*t-1)*ci^(s+1)*(n-m-2*t) - s*si^(n-m*2*t+1)*ci^(s-1));
                end
            end
        end
        
        wr = wr + -mu*Re^n/r^(n+1)*wrts;
        
    end
end

wr = 1/(mm*a*b)*wr;

% Osculation constraint: V = Rdot
wt = 0;

% Normal component
wn = mm*(1 + e*cos(f))^2/eta^3;

w = [wr, wt, wn];

% Test

% J2
Jcoeff = -Cnm(2:end,1);
w12 = 1/(mm*a*b)*Jcoeff(2)*mu*Re^2/r^3*3*sv*si*ci;


end

