function sw = BLanalyticalCoreyFoam(mug, muw)
%BLANALYTICAL Analytical solution of Buckley-Leverett equation
%muw = 10e-3;
%muo = 10e-3;
u = 1e-6/(pi()*0.038^2/4)/60;
phi = 0.2;
k = 2e-12; %[m^2]
sw_end = .5;
fmmob = 1.01889e+05;
fmdry = 1.20000e-01;
epdry = 5.0000e+03;
swc = 0.07;
sw_in = swc;
sgr = 0.0;
krg0 = 1;
ng = 2;
krw0 = 1;
nw = 2;
sws = @(sw)((sw>swc).*(sw-swc)/(1-sgr-swc));
kr = @(sw)(krg0*(1-sws(sw)).^ng);
% fm = @(sw)((sw>fmdry).*(1+fmmob*(0.5+atan(epdry.*(sw-fmdry))/pi()))+(sw<=fmdry));
fm = @(sw)(1+fmmob*(0.5+atan(epdry.*(sw-fmdry))/pi()));
krg = @(sw)(kr(sw)./fm(sw));
krw = @(sw)(krw0*sws(sw).^nw);
dkrwdsw = @(sw)(nw*krw0*(1/(1-sgr-swc))*sws(sw).^(nw-1));
dkrdsw = @(sw)((krg0*ng*(1-sws(sw)).^(ng-1))/(-swc-sgr+1));
% fm = @(sw)((sw>fmdry).*(1+fmmob*(0.5+atan(epdry.*(sw-fmdry))/pi()))+(sw<=fmdry));
dfmdsw = @(sw)(((epdry*fmmob)./(pi*(epdry^2*(sw-fmdry).^2+1))));
% dfmdsw = @(sw)((epdry*fmmob)./(pi*(epdry^2*(sw-fmdry).^2+1)));
dkrgdsw = @(sw)((dkrdsw(sw).*fm(sw)-dfmdsw(sw).*kr(sw))./fm(sw).^2);
fw = @(sw)((krw(sw)/muw)./(krw(sw)/muw+krg(sw)/mug));
dfwdsw = @(sw)((dkrwdsw(sw)/muw.*(krw(sw)/muw+krg(sw)/mug)- ...
    (dkrwdsw(sw)/muw+dkrgdsw(sw)/mug).*krw(sw)/muw)./ ...
    (krg(sw)/mug+krw(sw)/muw).^2);
s = linspace(swc,1,1000);
figure(1)
subplot(2,2,1);
plot(s, krw(s), s, kr(s), s, krg(s));
xlabel('S_w'); ylabel('rel perms'); legend('liq', 'gas', 'foam');
F = @(sw)(dfwdsw(sw)-(fw(sw_end)-fw(sw))/(sw_end-sw));
sw_shock = fzero(F, [swc+eps,0.3]);
subplot(2,2,2);
plot(s, fw(s), [sw_shock sw_end], [fw(sw_shock) fw(sw_end)]);
xlabel('S_w'); ylabel('f_w');
% plot(s, fw(s), [sw_end sw_shock], [fw(sw_end) fw(sw_shock)]);
s1 = linspace(sw_in, sw_shock, 50);
xt_s1 = u/phi*dfwdsw(s1);
xt_s = u/phi*dfwdsw(s);
xt_shock = u/phi*dfwdsw(sw_shock);
subplot(2,2,3);
plot(xt_s, s, '--', ...
    [xt_s1 xt_shock xt_shock max(xt_s)], [s1 sw_shock sw_end sw_end]) 
xlabel('x/t [m/s]'); ylabel('S_w');
subplot(2,2,4);
plot([xt_s1 xt_shock xt_shock 10*xt_shock], [s1 sw_shock sw_end sw_end])
xlabel('x/t [m/s]'); ylabel('S_w');
sw_prf = [s1 sw_shock sw_end sw_end];
xt_prf = [xt_s1 xt_shock xt_shock max(xt_s)];
sw = @(xt)([interp1(xt_s1, s1, xt(xt<xt_shock)) sw_end*xt(xt>=xt_shock)]);
% Now I can use the data to calculate the pressure drop across my domain
L = 0.17; % [m]
t = eps:20:5*40*60; % [s] time
dp = zeros(1,length(t));
for i =1:length(t)
    x = 0:0.001:L;
    xt = x/t(i);
    dp(i) = trapz(x, u./(k*(krg(sw(xt))/mug+krw(sw(xt))/muw)));
    
end
figure(2);plot(t/60, dp/1e5, '.')
xlabel('time [min]'); ylabel('pressure drop [bar]');
end

