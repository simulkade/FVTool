function [m_CO2, y_CO2] = DuanSun(T, p, m_salt)
% This function solves the Duan-Sun model to calculate the solubility of
% CO2 in water and in the brine and the mole fraction of CO2 in the gas phase
% T is the temperature in [K]
% p is the pressure in [Pa]
% m is the salt molality in [mol/kg]
% Copyright (c) 2012-2016 Ali Akbar Eftekhari
% See the license file
% R = 0.08314467; % [bar.L/(mol.K)]
P = p/1e5; % convert Pa to bar
% Constants (Table 2):
% constants for the calculation of chemical potential of CO2 in the liquid phase
c_mu_CO2_L = [28.9447706, -0.0354581768, -4770.67077, 1.02782768e-5, ...
		  33.8126098, 9.04037140e-3, -1.14934031e-3, -0.307405726, ...
		  -0.0907301486, 9.32713393e-4, 0.0];
% Constants for the calculation of second and third order interaction parameters
c_lambda_CO2_Na = [-0.411370585, 6.07632013e-4, 97.5347708, 0.0, 0.0, ...
				   0.0, 0.0, -0.0237622469, 0.0170656236, 0.0, 1.41335834e-5];
c_zeta_CO2_NaCl = [3.36389723e-4, -1.98298980e-5, 0.0, 0.0, 0.0, ...
				    0.0, 0.0, 2.12220830e-3, -5.24873303e-3, 0.0, 0.0];

% using the above constants (Eq. 7)
par = @(c)(c(1)+c(2)*T+c(3)/T+c(4)*T^2+c(5)/(630-T)+c(6)*P+ ...
		   c(7)*P*log(T)+c(8)*P/T+c(9)*P/(630-T)+ ...
		   c(10)*P^2/(630-T)^2+c(11)*T*log(P));

mu_CO2_L = par(c_mu_CO2_L);
lambda_CO2_Na = par(c_lambda_CO2_Na);
zeta_CO2_NaCl = par(c_zeta_CO2_NaCl);

% activity coefficient of CO2 (Eq. 5)
m_Na = m_salt;
m_Cl = m_salt;
log_gama_CO2 = 2*lambda_CO2_Na*m_Na+zeta_CO2_NaCl*m_Na*m_Cl;


% calculation of the fugacity of CO2 using the equation of Duan
a = [8.99288497e-2, -4.94783127e-1, 4.77922245e-2, 1.03808883e-2, ...
     -2.82516861e-2, 9.49887563e-2, 5.20600880e-4, -2.93540971e-4, ...
     -1.77265112e-3, -2.51101973e-5, 8.93353441e-5, 7.88998563e-5, ...
     -1.66727022e-2, 1.39800000,     2.96000000e-2];
% CO2 critical properties
Tc = 304.25; % [K]
Pc = 7.39e6; % [Pa]
Vc = 8.314467*Tc/Pc;
Tr = T/Tc;
Pr = P*1.0e5/Pc;
Z = @(Vr)(1+(a(1)+a(2)/Tr^2+a(3)/Tr^3)/Vr + ...
          (a(4)+a(5)/Tr^2+a(6)/Tr^3)/Vr^2 + (a(7)+a(8)/Tr^2+a(9)/Tr^3)/Vr^4 ...
          + (a(10)+a(11)/Tr^2+a(12)/Tr^3)/Vr^5 + ...
          a(13)/(Tr^3*Vr^2)*(a(14)+a(15)/Vr^2)*exp(-a(15)/Vr^2));
Vr_guess = (8.314*T*0.9/(P*1.0e5))/Vc;
F = @(Vr)(Z(Vr)/Vr-Pr/Tr);
Vr = fzero(F, Vr_guess);

% fugacity coefficient of pure CO2
Z_num = Z(Vr);
log_phi_CO2 = Z_num - 1.0 - log(Z_num) + (a(1)+a(2)/Tr^2+a(3)/Tr^3)/Vr + ...
			(a(4)+a(5)/Tr^2+a(6)/Tr^3)/(2.0*Vr^2) + ...
			(a(7)+a(8)/Tr^2+a(9)/Tr^3)/(4.0*Vr^4) + ...
			(a(10)+a(11)/Tr^2+a(12)/Tr^3)/(5.0*Vr^5) + ...
			a(13)/(2.0*Tr^3*a(15))*(a(14)+1.0-(a(14)+1.0+a(15)/Vr^2)*exp(-a(15)/Vr^2));

% pure water vapor pressure
c = [-38.640844, 5.8948420, 59.876516, 26.654627, 10.637097];
Pcw = 220.85; % [bar]
Tcw = 647.29; % [K]
t = (T-Tcw)/Tcw;
P_water = Pcw*T/Tcw*(1+c(1)*(-t)^1.9+c(2)*t+c(3)*t^2+c(4)*t^3+c(5)*t^4);
y_CO2 = (P-P_water)/P;

log_coef = mu_CO2_L - log_phi_CO2 + log_gama_CO2;
m_CO2 = y_CO2*P/exp(log_coef);
