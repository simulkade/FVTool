% capillary pressure curve for drainage
% Written by Ali A. Eftekhari
function res=dpc_drain(sw, pce, swc, labda)
pc0=1.0e7;
res=zeros(size(sw));
for i=1:numel(sw)  
  sw0=swc(i)+(1-labda*log(pc0/pce(i))+sqrt((-1+labda*log(pc0/pce(i)))^2+...
      4*swc(i)/(1-swc(i))))/2*(1-swc(i));
  if sw(i)>sw0
    res(i)=-1.0/((1-swc(i))*labda)*pce(i)*((sw(i)-swc(i))/(1-swc(i)))^(-1.0/labda-1);
  elseif 0.0<=sw(i) && sw(i)<=sw0
    res(i)=-1.0/((1-swc(i))*labda)*pce(i)*((sw0-swc(i))/(1-swc(i)))^(-1.0/labda-1);
  else
    res(i)=0.0;
  end
end