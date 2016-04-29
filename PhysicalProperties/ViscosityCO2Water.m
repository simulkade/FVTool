function [mu, muw]= ViscosityCO2Water(T, p, x_CO2)
% This function calculates the viscosity of the pure water, brine, and
% water-CO2, and brine-CO2 systems
% p is in Pa, T in K, x_CO2 is the mole fraction of CO2, mu in Pa.s
% Copyright (c) 2012-2016 Ali Akbar Eftekhari
% See the license file
Mwater = 0.0180153; % kg/mol
% MCO2 = 0.04401; % kg/mol

P = p/1e6; % convert Pa to MPa
m_CO2 = DuanSun(T, p, 0); % CO2 equilibrium molality (mol CO2/kg water)
x_CO2_eq = m_CO2/(m_CO2+1/Mwater);
muw = IAPWS_IF97('mu_pT', P, T); % pure water viscosity in Pa.s
vis_ratio = 1+(-4.069e-3*(T-273.15)+0.2531)*x_CO2/x_CO2_eq;
mu = muw*vis_ratio;
