function P = P_given_TVN(x)
% Downloaded from the personal page of Professor Sigurd Skogestad
% http://www.nt.ntnu.no/users/skoge/diplom/prosjekt03/wold/MatLab/
A = co2_sw(x);
P = -A.g(2);
