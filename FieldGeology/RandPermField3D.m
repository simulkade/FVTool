function [perm_field, poros_field, perm, poros] = RandPermField3D(k_avrg, phi_avrg, V_dp, m, n, w, M, N, W)
% generates permeability field on a m x n x w matrix, using the Dykstra-Parsons coefficient.
% the porousity field is generated based on the permeability value using the Carman-Kozeny equation
% we first build a log normal permeability and porousity field over m x n x w
% and then interpolate it over M x N x W
% Note that M > m, ...
s = -log(1-V_dp);
mu = log(k_avrg)-s^2/2;
% perm = exp(norminv(rand(m,n), mu, s));
% perm = exp(normrnd(mu, s, m, n));
perm = lognrnd(mu, s, m, n, w);
C = perm/k_avrg*phi_avrg^3/(1-phi_avrg)^2;
poros = zeros(m,n,w);
for i = 1:m
    for j = 1:n
        for k = 1:w
            r = roots([1 -C(i,j, k) 2*C(i,j,k) -C(i,j,k)]);
            poros(i,j,k) = min(r(imag(r)==0));
        end
    end
end

% interpolation
dx = 1/(m-1); dy = 1/(n-1); dz = 1/(w-1);
[X, Y, Z] = ndgrid(0:dx:1, 0:dy:1, 0:dz:1);
dxq = 1/(M-1); dyq = 1/(N-1); dzq = 1/(W-1);
[Xq, Yq, Zq] = ndgrid(0:dxq:1, 0:dyq:1, 0:dzq:1);

% P = [2 1 3];
% X = permute(X, P);
% Y = permute(Y, P);
% Z = permute(Z, P);
% perm = permute(perm, P);
% poros = permute(poros, P);
   
perm_field = interpn(X,Y,Z,perm,Xq,Yq,Zq, 'cubic');
perm_field(perm_field<=0) = interpn(X,Y,Z,perm,Xq(perm_field<=0),Yq(perm_field<=0),Zq(perm_field<=0), 'linear');
poros_field = interpn(X,Y,Z,poros,Xq,Yq,Zq, 'cubic');
poros_field(poros_field<=0) = interpn(X,Y,Z,poros,Xq(poros_field<=0),Yq(poros_field<=0),Zq(poros_field<=0), 'linear');
