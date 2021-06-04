% imbibition capillary pressure curve
% Written by Ali A. Eftekhari
function res=pc_imb(sw, pce, swc, sor, teta, labda, b)
  pc1=pc_drain(sw, pce, swc, labda);
  pc2=pc_drain(1-sw, pce, sor, labda);
  res=(0.5*(1+cos(teta))).^b.*pc1-(0.5*(1-cos(teta))).^b.*pc2;
end