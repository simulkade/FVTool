% imbibition capillary pressure curve
function res=dpc_imb(sw, pce, swc, sor, teta, labda, b)
  dpc1=dpc_drain(sw, pce, swc, labda);
  dpc2=dpc_drain(1-sw, pce, sor, labda);
  res=(0.5*(1+cos(teta))).^b.*dpc1+(0.5*(1-cos(teta))).^b.*dpc2;
end