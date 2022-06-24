function res = krw(sw, krw0, sor, swc, nw)
%KRW Corey-type water relative permeability
sws=sw_normalized(sw, swc, sor);
res=((sw<=1-sor).*krw0.*sws.^nw+(sw>1-sor).*(-(1-krw0)./sor.*(1.0-sw)+1.0));
end

