% relperm and capillary pressure curves
% to be used with FVTool package
% Written by Ali A. Eftekhari
% This is what this function does:
% res=zeros(size(sw));
% for i=1:numel(sw)  
%   if swc(i)<=sw(i) && sw(i)<=1-sor(i)
%     res(i)=kro0(i)*((1-sw(i)-sor(i))/(1-sor(i)-swc(i)))^no(i);
%   elseif 0.0<sw(i) && sw(i)<swc(i)
%     res(i)=1+(kro0(i)-1)/swc(i)*sw(i);
%   elseif sw(i)>1-sor(i)
%     res(i)=0.0;
%   elseif sw(i)<=0.0
%     res(i)=1.0;
%   end
% end
%
% However, it is slow in this form. Therefore, I have rewritten it in
% vectorized from. This one is kept here for reference. You can accelerate
% it by @parfor loop, but the vectorized form is substantially faster
function res=kro_loop(sw, kro0, sor, swc, no)
res=zeros(size(sw));
for i=1:numel(sw)  
  if swc(i)<=sw(i) && sw(i)<=1-sor(i)
    res(i)=kro0(i)*((1-sw(i)-sor(i))/(1-sor(i)-swc(i)))^no(i);
  elseif 0.0<sw(i) && sw(i)<swc(i)
    res(i)=1+(kro0(i)-1)/swc(i)*sw(i);
  elseif sw(i)>1-sor(i)
    res(i)=0.0;
  elseif sw(i)<=0.0
    res(i)=1.0;
  end
end