function h=h_given_TVN(x)
% Calculates spesific entalphy from given temperature, volume and
% mole number. Input x=[T;V;N]
% Downloaded from the personal page of Professor Sigurd Skogestad
% http://www.nt.ntnu.no/users/skoge/diplom/prosjekt03/wold/MatLab/
A = co2_sw(x);
h = [x(1),x(3)]*[-A.g(1);A.g(3)]/x(3);
