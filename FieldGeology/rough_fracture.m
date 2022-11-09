function d=rough_fracture(Nx,Ny, d_avg, sd, clx, cly, sd_cut_coef)
% All credits go to http://www.mysimlabs.com/surface_generation.html
% d_avg [m] average aperture between two rough surfaces of the fracture
% sd [m] standard deviation
Lx=1.0; % domain length
x=linspace(-Lx/2.0,Lx/2.0,Nx);
y=linspace(-Lx/2.0,Lx/2.0,Ny);
[X,Y]=ndgrid(x,y);
Z=sd*randn(Nx,Ny); % normal distribution
F = exp(-(X.^2/(clx*clx/2.0)+Y.^2/(cly*cly/2.0)));
% Gaussian filter
f =2.0/sqrt(pi)*Lx/sqrt(Nx*Ny)/sqrt(clx)/sqrt(cly).*ifft2(fft2(Z).*fft2(F)); % another filter
d=d_avg+real(f); % aperture field
d(d<0)=0.0;
d(d>sd_cut_coef*sd) = sd_cut_coef*sd;

