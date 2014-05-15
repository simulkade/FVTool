function c_CO2 = CO2GasDensity(T, p)
% This function uses the Span and Wagner equation of state for the pure CO2
% and calculates the molar density of CO2 in mol/m^3 as a function of
% temperature (K) and pressure (Pa)
% I have not written the Span-Wagner EOS myself.

c_CO2 = 1/V_given_TPN([T;p;1],'v'); %mol/m^3