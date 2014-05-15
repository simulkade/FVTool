function V = V_given_TPN(x, phase)
% Calculates volume for given temperature, pressure and mole number
% Usage :  V = V_given_TPN(x,'phase')
% Inputs:   x       = [T;P;N]
%           phase   = l (liquid) or
%                   = v (vapour)
% Downloaded from the personal page of Professor Sigurd Skogestad
% http://www.nt.ntnu.no/users/skoge/diplom/prosjekt03/wold/MatLab/
R = 8.314;
T = x(1);
P = x(2);
N = x(3);
if phase=='v'
  V = N*R*T/P;
elseif phase=='l'
  V = N*R*T/(10*P);
end
if P>80e5
    V = N*R*T/(10*P);
end

Pt = P+1e9;
dx = 0;
while abs(P-Pt)/P>1e-6
  V=V+dx;
  A=co2_sw([T;V;N]);
  Pt=-A.g(2);
  dx = -inv(A.H(2,2))*(P-Pt);
end
    
