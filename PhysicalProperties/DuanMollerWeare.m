function Vr = DuanMollerWeare(Tr,Pr)
% This is the EOS of Duan-Moller-Weare for the CO2-Methane-water system
% This version calculates the density of the CO2-water and the
% CO2-water-NaCl system
% Written by Ali A. Eftekhari at TU Delft
% See the license file

% constants
a_CH4 = [8.72553928e-2, -7.52599476e-1, 3.75419887e-1, 1.07291342e-2, ...
    5.49626360e-3, -1.84772802e-2, 3.18993183e-4, 2.11079375e-4, ...
    2.01682801e-5, -1.65606189e-5, 1.19614546e-4, -1.08087289e-4];
alfa_CH4 = 4.48262295e-02;  beta_CH4 = 7.5397e-01; gama_CH4 = 7.7167e-02;

a_CO2 = [8.99288497e-02, -4.94783127e-01, 4.77922245e-02, 1.03808883e-02, ...
    -2.82516861e-02, 9.49887563e-02, 5.20600880e-04, -2.93540971e-04, ...
    -1.77265112e-03, -2.51101973e-05, 8.93353441e-05, 7.88998563e-05];
alfa_CO2 = -1.66727022e-02;  beta_CO2 = 1.398; gama_CO2 = 2.96e-02;

a_water = [8.64449220e-02, -3.96918955e-01, -5.73334886e-02, -2.93893000e-04, ...
    -4.15775512e-03, 1.99496791e-02, 1.18901426e-04, 1.55212063e-04, ...
    -1.06855859e-04, -4.93197687e-06, -2.73739155e-06, 2.65571238e-06];
alfa_CO2 = 8.96079018e-03;  beta_CO2 = 4.02; gama_CO2 = 2.57e-02;

B = @(a,Tr)(a(1)+a(2)/Tr^2+a(3)/Tr^3);
C = @(a,Tr)(a(4)+a(5)/Tr^2+a(6)/Tr^3);
D = @(a,Tr)(a(7)+a(8)/Tr^2+a(9)/Tr^3);
E = @(a,Tr)(a(10)+a(11)/Tr^2+a(12)/Tr^3);
F = @(alfa,Tr)(alfa/Tr^3);
Z = @(Vr,a, alfa_c, beta_c, gama_c, Tr)(1+B(a,Tr)/Vr+C(a,Tr)/Vr^2+ ...
    D(a,Tr)/Vr^4+E(a,Tr)/Vr^5+ ...
    F(alfa,Tr)/Vr^2*(beta_c+gama_c/Vr^2)*exp(-gama_c/Vr^2));
