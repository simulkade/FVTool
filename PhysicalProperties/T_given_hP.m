function T = T_given_hP(x,h,P)
% Downloaded from the personal page of Professor Sigurd Skogestad
% http://www.nt.ntnu.no/users/skoge/diplom/prosjekt03/wold/MatLab/

htmp = 10*h;
Ptmp = 10*P;
while norm([(h-htmp)/h;(P-Ptmp)/P])>1e-8
  A = co2_sw(x);
  Ptmp = -A.g(2);
  htmp = [x(1),x(3)]*[-A.g(1);A.g(3)]/x(3);
  x(1) = x(1)-100*(h-htmp)*x(3)/(-(A.g(1)^2+x(1)*A.H(1,1))+x(3)*A.H(3,1));
  x(2) = x(2)-(P-Ptmp)/(A.H(2,2));
end
T = x(1);
