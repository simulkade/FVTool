function res = dkrwdsw(sw, krw0, sor, swc, nw)
%KRW Corey-type water relative permeability
sws=sw_normalized(sw, swc, sor);
res=((sw<=1-sor).*nw.*krw0.*(1./(1-sor-swc)).*sws.^(nw-1)+(sw>1-sor).*((1-krw0)./sor));
end

