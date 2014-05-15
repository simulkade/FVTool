function molar_dens = JacobsenStewart_N2(p_pa, T)
%JACOBSENSTEWART_N2 calculates the molar density of N2 as a function of
%pressure (Pa) and temperature (K) using the equation proposed by Jacobsen
%and Stewart. The equation is valid from 63-2000 K with pressures up to 10000
%bar
%   Detailed explanation goes here
% molar_dens: [mol/m^3] gas phase molar density

p = p_pa/101325; % convert Pa to atm

% constants
N = [0.136097243686983e-2
     0.107028500555504
    -0.243926251659031e1
     0.341240789637052e2
    -0.422956791527436e4
     0.105277159433708e-3
    -0.111355267180312e-1
     0.142748464727047e-3
     0.179621096187021e5
     0.751267113751007e-7
     0.231737284741220e-2
    -0.509008258448481
     0.488523311385766e-4
    -0.112001704676209e-2
    -0.678366343173806
     0.742796115735318e-4
    -0.110400310752087e-5
     0.341309483327025e-3
    -0.166216790652177e-5
    -0.164616585853003e5
    -0.119724198386865e6
    -0.948085610750225e2
     0.554879602331972e5
    -0.174677685666810
    -0.256709963280944e1
    -0.404528219006087e-3
    -0.257279422571519
    -0.121204517440575e-6
     0.104690038752288e-4
    -0.529499586313775e-9
    -0.774723053052639e-8
     0.610368224362452e-7];
 R = 0.0820562; % liter.atm/(mol.K)
 labda = 0.0056;
 
 f = @(rho)(-p+rho*R*T ...
        + rho^2*(N(1)*T+N(2)*T^0.5+N(3)+N(4)/T+N(5)/T^2)+ ...
        + rho^3*(N(6)*T+N(7)+N(8)/T+N(9)/T^2) ...
        + rho^4*(N(10)*T+N(11)+N(12)/T) + rho^5*N(13) ...
        + rho^6*(N(14)/T+N(15)/T^2) + rho^7*N(16)/T ...
        + rho^8*(N(17)/T+N(18)/T^2) + rho^9*N(19)/T^2 ...
        + rho^3*(N(20)/T^2+N(21)/T^3)*exp(-labda*rho^2) ...
        + rho^5*(N(22)/T^2+N(23)/T^4)*exp(-labda*rho^2) ...
        + rho^7*(N(24)/T^2+N(25)/T^3)*exp(labda*rho^2) ...
        + rho^9*(N(26)/T^2+N(27)/T^4)*exp(-labda*rho^2) ...
        + rho^11*(N(28)/T^2+N(29)/T^3)*exp(-labda*rho^2) ...
        + rho^13*(N(30)/T^2+N(31)/T^3+N(32)/T^4)*exp(-labda*rho^2));

rho0 = p/(R*T); % mol/liter initial estimate
molar_dens = fzero(f, rho0)*1000; % mol/m^3
