function [xt_prf, sw_prf, rec_fact] = BLanalyticalCorey(muw, muo)
%BLANALYTICAL Analytical solution of Buckley-Leverett equation
% Copyright (c) 2012-2016 Ali Akbar Eftekhari
% See the license file
%muw = 10e-3;
%muo = 10e-3;
pv_inj=5; % injected pore volume
L=1.0; % m core length
u = 1e-4;
phi = 0.2;
sw_in = 1;
krw0 = 1.0;
kro0 = 0.76;
nw = 2.4;
no = 2.0;
sor=0.12;
swc=0.09;
sw_end = swc;
sws=@(sw)((sw>swc).*(sw<1-sor).*(sw-swc)/(1-sor-swc)+(sw>=1-sor).*ones(size(sw)));
kro=@(sw)((sw>=swc).*kro0.*(1-sws(sw)).^no+(sw<swc).*(1+(kro0-1)/swc*sw));
krw=@(sw)((sw<=1-sor).*krw0.*sws(sw).^nw+(sw>1-sor).*(-(1-krw0)/sor.*(1.0-sw)+1.0));
dkrwdsw=@(sw)((sw<=1-sor).*nw.*krw0.*(1/(1-sor-swc)).*sws(sw).^(nw-1)+(sw>1-sor)*((1-krw0)/sor));
dkrodsw=@(sw)((sw>=swc).*(-kro0*no*(1-sws(sw)).^(no-1))/(-swc-sor+1)+(sw<swc).*((kro0-1)/swc));
fw = @(sw)((krw(sw)/muw)./(krw(sw)/muw+kro(sw)/muo));
dfwdsw = @(sw)((dkrwdsw(sw)/muw.*(krw(sw)/muw+kro(sw)/muo)- ...
    (dkrwdsw(sw)/muw+dkrodsw(sw)/muo).*krw(sw)/muw)./ ...
    (krw(sw)/muw+kro(sw)/muo).^2);
s = linspace(0,1,100);
figure(1)
subplot(2,2,1);
plot(s, krw(s), s, kro(s));
F = @(sw)(dfwdsw(sw)-(fw(sw)-fw(sw_end))/(sw-sw_end));
sw_shock = fzero(F, [swc+eps,1-sor-eps]);
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
xt_prf = [xt_s1 xt_shock xt_shock+eps 2*xt_shock];

% make data better!
i=1;
while(true)
  if (i+1)==length(xt_prf)
    break
  elseif xt_prf(i)>=xt_prf(i+1)
    xt_prf(i+1)=[];
    sw_prf(i+1)=[];
  else
    i=i+1;
  end
end
% calculate the recovery factor
x = linspace(0,L,1000);
xt_int=[xt_prf L/eps];
sw_int=[sw_prf sw_end];
t_inj=pv_inj*phi*L/u;
t = linspace(eps(),t_inj, 100); % [s] time
rec_fact = zeros(length(t),1);
rec_fact(1)=0.0;
for k = 2:length(t)
    xt = x/t(k);
    swt=interp1(xt_int, sw_int, xt);
    rec_fact(k) = 1-trapz(x, 1-swt);
end
subplot(2,2,4);
plot(t/(phi*L/u), rec_fact); xlabel('PV'); ylabel('recovery factor');
end
