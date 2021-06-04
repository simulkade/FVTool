function [xt_prf, sw_prf] = BLanalytical(muw, muo)
%BLANALYTICAL Analytical solution of Buckley-Leverett equation
% Written by Ali A. Eftekhari
% See the license file
%muw = 10e-3;
%muo = 10e-3;
u = 1e-3;
phi = 0.2;
sw_in = 1;
sw_end = 0;
krw = @(sw)(sw.^4);
dkrwdsw = @(sw)(4*sw.^3);
kro = @(sw)((1-sw.^2).*(1-sw).^2);
dkrodsw = @(sw)(-2*sw.*(1-sw).^2-2*(1-sw).*(1-sw.^2));
fw = @(sw)((krw(sw)/muw)./(krw(sw)/muw+kro(sw)/muo));
dfwdsw = @(sw)((dkrwdsw(sw)/muw.*(krw(sw)/muw+kro(sw)/muo)- ...
    (dkrwdsw(sw)/muw+dkrodsw(sw)/muo).*krw(sw)/muw)./ ...
    (krw(sw)/muw+kro(sw)/muo).^2);
s = linspace(0,1,100);
figure(1)
subplot(2,2,1);
plot(s, krw(s), s, kro(s));
F = @(sw)(dfwdsw(sw)-(fw(sw)-fw(sw_end))/(sw-sw_end));
sw_shock = fzero(F, [eps,1]);
subplot(2,2,2);
plot(s, fw(s), [sw_end sw_shock], [fw(sw_end) fw(sw_shock)]);
% plot(s, fw(s), [sw_end sw_shock], [fw(sw_end) fw(sw_shock)]);
s1 = linspace(sw_in, sw_shock, 50);
xt_s1 = u/phi*dfwdsw(s1);
xt_s = u/phi*dfwdsw(s);
xt_shock = u/phi*dfwdsw(sw_shock);
subplot(2,2,3);
plot(xt_s, s, '--', ...
    [xt_s1 xt_shock xt_shock max(xt_s)], [s1 sw_shock sw_end sw_end])
sw_prf = [s1 sw_shock sw_end sw_end];
xt_prf = [xt_s1 xt_shock xt_shock max(xt_s)];
end
