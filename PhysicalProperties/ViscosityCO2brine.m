function [mu, muw, mu_brine]= ViscosityCO2brine(T, p, x_CO2, m_salt)
% This function calculates the viscosity of the pure water, brine, and
% water-CO2, and brine-CO2 systems
% p is in Pa, T in K, x_CO2 is the mole fraction of CO2, mu in Pa.s
% Written by Ali A. Eftekhari at TU Delft
% See the license file
Mwater = 0.0180153; % kg/mol
% MCO2 = 0.04401; % kg/mol

P = p/1e6; % convert Pa to MPa
m_CO2 = DuanSun(T, p, 0); % CO2 equilibrium molality (mol CO2/kg water)
x_CO2_eq = m_CO2/(m_CO2+1/Mwater);
muw = IAPWS_IF97('mu_pT', P, T); % pure water viscosity in Pa.s

a0 = -0.21319213; a1 = 0.13651589e-2;  a2 = -0.12191756e-5;
b0 = 0.69161945e-1; b1 = -0.27292263e-3; b2 = 0.20852448e-6;
c0 = -0.25988855e-2;  c1 = 0.77989227e-5;

A = a0 + a1*T+a2*T*T;
B = b0 + b1*T+b2*T*T;
C = c0 + c1*T;

mu_brine = muw*exp(A*m_salt+B*m_salt^2+C*m_salt^3);

vis_ratio = 1+(-4.069e-3*(T-273.15)+0.2531)*x_CO2/x_CO2_eq;
mu = mu_brine*vis_ratio;
