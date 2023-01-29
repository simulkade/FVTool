% capillary pressure curve for drainage
% Written by Ali A. Eftekhari
function res=pc_drain(sw, pce, swc, labda, pc_max)
sw0=swc+(1-labda*log(pc_max/pce)+sqrt((-1+labda*log(pc_max/pce))^2+...
      4*swc/(1-swc)))/2*(1-swc);
pcs=pce*((sw0-swc)/(1-swc)).^(-1.0/labda);
res = (sw>sw0).*(pce*((sw-swc+eps)/(1-swc)).^(-1.0/labda))+...
(0.0<=sw).*(sw<=sw0).*(exp((log(pcs)-log(pc_max))/sw0*(sw-sw0)+log(pcs)))+...
(sw<0)*pc_max;
end