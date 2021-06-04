% the value of water saturation at pc=0
% Written by Ali A. Eftekhari
function res=sw_zero_pc_imb(swc, sor, teta, labda, b)
    if teta~=pi
        res=((1-swc)^2+swc*(1-sor)*((1-cos(teta))/(1+cos(teta)))^(b*labda))...
            /((1-swc)+(1-sor)*((1-cos(teta))/(1+cos(teta)))^(b*labda));
    else
        res=swc;
    end
end