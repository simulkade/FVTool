function [S] = co2_sw(x)
% co2_sw.m
% Downloaded from the personal page of Professor Sigurd Skogestad
% http://www.nt.ntnu.no/users/skoge/diplom/prosjekt03/wold/MatLab/
[S] = Helmholtz_anonymous_23829804(x);


function [S] = Helmholtz_anonymous_23829804(x)

[S] = StandardState_anonymous_23817360(x);


function [S] = StandardState_anonymous_23817360(x)

[S]      = EquationOfState_anonymous_23031408(x);
i        = [3];
[S_1]    = MuT_cp_dippr_23809380(x(1));
S.g(1)   = S.g(1) + S_1.dmudT'*x(i);
S.g(i)   = S.g(i) + S_1.mu;
S.H(1,1) = S.H(1,1) + S_1.d2mudTdT'*x(i);
S.H(i,1) = S.H(i,1) + S_1.dmudT;
S.H(1,i) = S.H(i,1)';


function [S] = EquationOfState_anonymous_23031408(x)

[S]   = EquationOfState_anonymous_22908308(x);
[S_1] = ModTVN_ideal_idealgas_23028984(x);
S.g   = S.g + S_1.g;
S.H   = S.H + S_1.H;


function [S] = EquationOfState_anonymous_22908308(x)

[S] = ModTVN_sw_0_5_7_0_23818404(x);


function [S] = ModTVN_sw_0_5_7_0_23818404(x)

R               = 8.314511984;
T_c             = 304.1282;
rho_c           = 10624.90627;
tau             = T_c/x(1);
delta           = x(3)/x(2)/rho_c;
N               = x(3);
a_1             = [0.89875108;-2.1281985;-0.06819032;0.076355306;0.00022053253];
t_1             = [0.25;1.25;1.5;0.25;0.875];
d_1             = [1.0;1.0;1.0;3.0;7.0];
a_2             = [0.41541823;0.71335657;0.00030354234;-0.36643143;-0.0014407781;-0.089166707;-0.023699887];
t_2             = [2.375;2.0;2.125;3.5;6.5;4.75;12.5];
d_2             = [1.0;2.0;5.0;1.0;1.0;4.0;2.0];
p_2             = [1.0;1.0;1.0;2.0;2.0;2.0;3.0];
u_2             = d_2 - p_2.*delta.^p_2;
v_2             = exp(-delta.^p_2);
phir            = a_1'*(delta.^d_1.*tau.^t_1) + a_2'*(v_2.*delta.^d_2.*tau.^t_2);
phi_taur        = a_1'*(delta.^d_1.*t_1.*tau.^(t_1 - 1)) + a_2'*(v_2.*delta.^d_2.*t_2.*tau.^(t_2 - 1));
phi_deltar      = a_1'*(d_1.*delta.^(d_1 - 1).*tau.^t_1) + a_2'*(v_2.*delta.^(d_2 - 1).*u_2.*tau.^t_2);
phi_tautaur     = a_1'*(delta.^d_1.*t_1.*(t_1 - 1).*tau.^(t_1 - 2)) + a_2'*(v_2.*delta.^d_2.*t_2.*(t_2 - 1).*tau.^(t_2 - 2));
phi_deltadeltar = a_1'*(d_1.*(d_1 - 1).*delta.^(d_1 - 2).*tau.^t_1) + a_2'*(v_2.*delta.^(d_2 - 2).*(u_2.*(u_2 - 1) - p_2.^2.*delta.^p_2).*tau.^t_2);
phi_deltataur   = a_1'*(d_1.*delta.^(d_1 - 1).*t_1.*tau.^(t_1 - 1)) + a_2'*(v_2.*delta.^(d_2 - 1).*t_2.*u_2.*tau.^(t_2 - 1));
g_1             = N*R*(phir - tau*phi_taur);
g_2             = -R*T_c*rho_c*delta^2*phi_deltar/tau;
g_i             = R*T_c*(phir + delta*phi_deltar)/tau;
H_11            = N*R*tau^3*phi_tautaur/T_c;
H_21            = -R*delta^2*rho_c*(phi_deltar - tau*phi_deltataur);
H_i1            = R*(phir - tau*phi_taur + delta*(phi_deltar - tau*phi_deltataur));
H_22            = R*T_c*rho_c^2*delta^3*(delta*phi_deltadeltar + 2*phi_deltar)/(tau*N);
H_i2            = -R*T_c*rho_c*delta^2*(delta*phi_deltadeltar + 2*phi_deltar)/(tau*N);
H_ii            = R*T_c*delta*(delta*phi_deltadeltar + 2*phi_deltar)/(tau*N);
S.g             = [g_1;g_2;g_i];
S.H             = [H_11,H_21,H_i1;H_21,H_22,H_i2;H_i1,H_i2,H_ii];


function [S] = ModTVN_ideal_idealgas_23028984(x)

R     = 8.314511984;
pcirc = [101325.0];
xcirc = [1.0];
T     = x(1);
V     = x(2);
i     = [3];
n     = x(i);
NR    = sum(n)*R;
NRT   = NR*T;
e     = [1];
p     = R*log(R*T*(n./pcirc)/V);
g_1   = n'*p;
g_2   = -NRT/V;
g_i   = T*p;
H_11  = NR/T;
H_21  = -NR/V;
H_i1  = p + R*e;
H_22  = NRT/V^2;
H_i2  = -R*(T/V)*e;
H_ii  = R*T*diag(e./n);
S.g   = [g_1;g_2;g_i];
S.H   = [H_11,H_21',H_i1';H_21,H_22,H_i2';H_i1,H_i2,H_ii];


function [S] = MuT_cp_dippr_23809380(T)

t_0        = [298.15];
t_min      = [50.0];
t_max      = [5000.0];
C          = [29.37,34.54,1428.0,26.4,588.0];
c_1        = C(:,1);
c_2        = C(:,2);
c_3        = C(:,3);
c_4        = C(:,4);
c_5        = C(:,5);
c_p        = c_1 + c_2.*c_3.^2./(T^2*sinh(c_3/T).^2) + c_4.*c_5.^2./(T^2*cosh(c_5/T).^2);
h          = c_1*T + c_2.*c_3.*cosh(c_3/T)./sinh(c_3/T) - c_4.*c_5.*sinh(c_5/T)./cosh(c_5/T) - c_1.*t_0 - c_2.*c_3.*cosh(c_3./t_0)./sinh(c_3./t_0) + c_4.*c_5.*sinh(c_5./t_0)./cosh(c_5./t_0);
s          = c_1*log(T) + c_4.*(log(cosh(c_5/T)) - c_5.*sinh(c_5/T)./cosh(c_5/T)/T) + c_2.*(c_3.*cosh(c_3/T)./sinh(c_3/T)/T - log(sinh(c_3/T))) - c_1.*log(t_0) - c_4.*(log(cosh(c_5./t_0)) - c_5.*sinh(c_5./t_0)./cosh(c_5./t_0)./t_0) - c_2.*(c_3.*cosh(c_3./t_0)./sinh(c_3./t_0)./t_0 - log(sinh(c_3./t_0)));
S.mu       = h - T*s;
S.dmudT    = -s;
S.d2mudTdT = -(c_p/T);
i          = [1];
[S_1]      = MuT_hs_h0_23796420(T);
S.mu(i)    = S.mu(i) + S_1.mu;
S.dmudT(i) = S.dmudT(i) + S_1.dmudT;


function [S] = MuT_hs_h0_23796420(T)

hcirc   = [-393510.0];
[S]     = MuT_hs_s0_23790264(T);
S.mu    = S.mu + hcirc;
S.dmudT = S.dmudT + [0];


function [S] = MuT_hs_s0_23790264(T)

scirc   = [213.677];
S.mu    = -(T*scirc);
S.dmudT = -scirc;

