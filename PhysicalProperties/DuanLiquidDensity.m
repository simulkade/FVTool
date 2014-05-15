function [rho, rhob] = DuanLiquidDensity(T, p, x_CO2, x_salt)
% This function uses the model of Duan, Hu, Li, and Mao to calculate the
% liquid density of the CO2-H2O and CO2-H2O-NaCl systems up to 647 K and
% 100 MPa
% T in K
% p in Pa
% x_CO2, x_salt: mole fractions of CO2 and NaCl in the liquid phase

% constants
% molecular mass (kg/mol):
Mwater = 0.0180153;
MCO2 = 0.04401;
MNaCl = 0.05844;


% Eq (7):
% P in Mpa, T in K, V1 is the molar volume of pure water (IAPWS97)
P = p/1e6; % convert p to MPa
V1 = IAPWS_IF97('v_pT', P, T)*Mwater; % m^3/mol
rhow = Mwater/V1; % water density (kg/m^3)
x2 = x_CO2;
A1j = [0.38384020e-3, -0.55953850, 0.30429268e3, -0.72044305e5, 0.63003388e7];
A2j = [-0.57709332e-5, 0.82764653e-2, -0.43813556e1, 0.10144907e4, -0.86777045e5];
A = @(Ai)(Ai(1)*T^2+Ai(2)*T+Ai(3)+Ai(4)/T+Ai(5)/T^2);
A1 = A(A1j); A2 = A(A2j);

if x_salt == 0
    V = V1*(1+(A1+A2*P)*x2); % m^3/mol
    M = x_CO2*MCO2 + (1-x_CO2)*Mwater; % kg/mol
    rho = M/V; % kg/m^3
    rhob = rhow;
else
    x_water = 1-x_salt-x_CO2; % water mole fraction
    y_water = x_water/(x_water+x_salt); y_salt = 1-y_water; % normalized brine mole fraction
    w_salt = y_salt*MNaCl/(y_salt*MNaCl+y_water*Mwater); % normalized salt mass fraction
    rhob = rhow + w_salt*1000*(0.668+0.44*w_salt+1e-6*(300*P-2400*P*w_salt + ...
        T*(80+3*T-3300*w_salt-13*P+47*P*w_salt))); % Batzle-Wang 1992
    Msw = y_water*Mwater+y_salt*MNaCl; % kg/mol normalized brine
    Vb = Msw/rhob; % m^3/mol brine molar volume
    V_phi3 = (Vb-y_water*V1)/y_salt; % apparent molar volume of NaCl in brine
    V_phi2 = V1*(1+A1+A2*P); % apparent molar volume of CO2 in ternary
    V = x_water*V1+x_CO2*V_phi2+x_salt*V_phi3; % molar volume m^3/mol
    M = x_CO2*MCO2 + x_water*Mwater+ x_salt*MNaCl;
    rho = M/V; % kg/m^3
end
    