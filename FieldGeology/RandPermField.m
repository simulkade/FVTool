function [perm_field, poros_field, perm, poros] = RandPermField(k_avrg, phi_avrg, V_dp, m, n, M, N)
% generates permeability field on a m x n matrix, using the Dykstra-Parsons coefficient.
% the porosity field is generated based on the permeability value using the Carman-Kozeny equation
% we first build a log normal permeability and porousity field over m x n
% and then interpolate it over M x N
s = -log(1-V_dp);
mu = log(k_avrg)-s^2/2;
% perm = exp(norminv(rand(m,n), mu, s));
% perm = exp(normrnd(mu, s, m, n));
perm = lognrnd(mu, s, m, n);
C = perm/k_avrg*phi_avrg^3/(1-phi_avrg)^2;
poros = zeros(m,n);
for i = 1:m
	for j = 1:n
		r = roots([1 -C(i,j) 2*C(i,j) -C(i,j)]);
		poros(i,j) = min(r(imag(r)==0));
	end
end

% interpolation
dx = 1/(m-1); dy = 1/(n-1);
[X, Y] = meshgrid(0:dx:1,0:dy:1);
dxq = 1/(M-1); dyq = 1/(N-1);
[Xq, Yq] = meshgrid(0:dxq:1, 0:dyq:1);

 log_perm_field = interp2(X,Y,log(perm'),Xq,Yq, 'cubic')';
 perm_field = exp(log_perm_field);
%  perm_field(perm_field<=0) = interp2(X,Y,perm', ...
%      Xq(perm_field<=0),Yq(perm_field<=0), 'linear');
 poros_field = interp2(X,Y,poros',Xq,Yq, 'cubic')';
 poros_field(poros_field<=0) = interp2(X,Y,poros', ...
     Xq(poros_field<=0),Yq(poros_field<=0), 'linear');
