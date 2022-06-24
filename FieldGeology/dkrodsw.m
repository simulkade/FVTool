% function res=dkrodsw(sw, kro0, sor, swc, no)
% derivative of Corey-type water relative permeability curve
function res=dkrodsw(sw, kro0, sor, swc, no)
sws=sw_normalized(sw, swc, sor);
res=((sw>=swc).*(-kro0.*no.*(1-sws).^(no-1))./(-swc-sor+1)+(sw<swc).*((kro0-1)./swc));
end