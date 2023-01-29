% capillary pressure curve for drainage
% Written by Ali A. Eftekhari
% TBD, probably not required in the new IMPES formulation
function res=dpc_drain(sw, pce, swc, labda, pc_max)
sw0=swc+(1-labda*log(pc_max/pce)+sqrt((-1+labda*log(pc_max/pce))^2+...
      4*swc/(1-swc)))/2*(1-swc);
res = (sw>sw0).*(-1.0/((1-swc)*labda)*pce*((sw-swc)/(1-swc)).^(-1.0/labda-1))+...
    (0.0<=sw).*(sw<=sw0).*(-1.0/((1-swc)*labda)*pce*((sw0-swc)/(1-swc)).^(-1.0/labda-1))+...
    (sw<0)*0.0;
end