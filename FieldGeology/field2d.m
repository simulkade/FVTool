function perm=field2d(Nx,Ny,k_avg,V_dp,clx,cly)
% Copyright (c) 2012-2016 Ali Akbar Eftekhari
% See the license file
Lx=1.0; % domain length
x=linspace(-Lx/2.0,Lx/2.0,Nx);
y=linspace(-Lx/2.0,Lx/2.0,Ny);
[X,Y]=ndgrid(x,y);
s=-log(1-V_dp); % standard deviation
mu=log(k_avg)-s*s/2.0; % mean for a log-random field
Z=s*randn(Nx,Ny); % normal distribution
F = exp(-(X.^2/(clx*clx/2.0)+Y.^2/(cly*cly/2.0)));
% Gaussian filter
f =2.0/sqrt(pi)*Lx/sqrt(Nx*Ny)/sqrt(clx)/sqrt(cly).*ifft2(fft2(Z).*fft2(F)); % another filter
perm=exp(mu+real(f)); % perm field
