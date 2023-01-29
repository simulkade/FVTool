% imbibition capillary pressure curve
% Written by Ali A. Eftekhari
function res=dpc_imb(sw, pce_w, pce_o, swc, sor, labda_w, labda_o, pc_max_w, pc_max_o)
  dpc1=dpc_drain(sw, pce_w, swc, labda_w, pc_max_w);
  dpc2=dpc_drain(1-sw, pce_o, sor, labda_o, pc_max_o);
  res=dpc1-dpc2;
end