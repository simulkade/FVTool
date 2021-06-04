function perm=field3d(Nx,Ny,Nz,k_avg,V_dp,clx,cly,clz)
% Written by Ali A. Eftekhari
% See the license file
% Does not work! Needs lots of improvements
Lx=1.0; % domain length
x=linspace(-Lx/2.0,Lx/2.0,Nx);
y=linspace(-Lx/2.0,Lx/2.0,Ny);
z=linspace(-Lx/2.0,Lx/2.0,Nz);
[X,Y,Z]=ndgrid(x,y,z);
s=-log(1-V_dp); % standard deviation
mu=log(k_avg)-s*s/2.0; % mean for a log-random field
Z=s*randn(Nx,Ny,Nz); % normal distribution
F = exp(-(X.^2/(clx*clx/2.0)+Y.^2/(cly*cly/2.0)+Z.^2/(clz*clz/2.0)));
% Gaussian filter
f =2.0/sqrt(pi)*Lx/(Nx*Ny*Nz)^(1/3)/(clx*cly*clz)^(1/3).*ifftn(fft2(Z).*fftn(F)); % another filter
perm=exp(mu+real(f)); % perm field
