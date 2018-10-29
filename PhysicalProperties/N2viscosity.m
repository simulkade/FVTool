function mu = N2viscosity(p_pa, T)
% this function calculates the viscosity of N2 as a function of temperature
% T [K] and pressure p [Pa] using the formulation of Lemmon and Jacobsen
% Copyright (c) 2012-2016 Ali Akbar Eftekhari
% See the license file
Tc = 126.192; % [K] critical temperature
rho_c = 11.1839; % mol/dm^3
% pc = 3.3958; % [MPa] critical pressure
M = 28.01348; % [g/mol]
eps_over_k = 98.94; % [K] eps/k Leonard-Jones energy
sigma = 0.3656; % [nm]
% zeta0 = 0.17; % [nm]
% GAMA = 0.055; % [-]
% qD = 0.40; % [nm]
% Tref = 252.384; % [K]

rho = JacobsenStewart_N2(p_pa, T)/1000; % mol/dm^3
delta = rho/rho_c;
tau = Tc/T;

ni = [0 1 2 3 4];
b = [0.431 -0.4623 0.08406 0.005341 -0.00331]; % fitted parameters
Tstar = T/eps_over_k;
RHO = exp(sum(b.*log(Tstar).^ni));
eta0 = 0.0266958*sqrt(M*T)/(sigma^2*RHO);

N = [10.72 0.03989 0.001208 -7.402 4.620];
t = [0.1 0.25 3.2 0.9 0.3];
d = [2 10 12 2 1];
l = [0 1 1 2 3];
gama = [0 1 1 1 1];

eta_r = sum(N.*tau.^t.*delta.^d.*exp(-gama.*delta.^l));



eta = eta0 + eta_r; % eta is viscosity in [micro Pa.s]

mu = eta*1e-6; % viscosity in Pa.s
