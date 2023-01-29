% imbibition capillary pressure curve
% Written by Ali A. Eftekhari
% Note that pc_max_o is specified as pc_min in the input json files for
% slightly more consistency witht the polynomial Pc
function res=pc_imb(sw, pce_w, pce_o, swc, sor, labda_w, labda_o, pc_max_w, pc_max_o)
  pc1=pc_drain(sw, pce_w, swc, labda_w, pc_max_w);
  pc2=pc_drain(1-sw, pce_o, sor, labda_o, pc_max_o);
  res=pc1-pc2;
end