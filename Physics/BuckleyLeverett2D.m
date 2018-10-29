function [sw, p] = BuckleyLeverett2D(MeshStructure, BCp, BCs, ...
			perm, mu_oil, mu_water, sw0, p0, dt, s, eps1)
% This function solves the continuity equation for an incompressible fluid
% flow in a porous medium, where Darcy's law is applicable.
% Copyright (c) 2012-2016 Ali Akbar Eftekhari
% See the license file 
% extract the mesh data
mn = MeshStructure.numberofcells;
m = mn(1); n = mn(2);



% start the solver here
% first solve the convection equation for sw
% initialize
% sw.Old = cellBoundary(MeshStructure, BCs, sw0);
FL = fluxLimiter('SUPERBEE');
sw.Old = sw0;
sw.value = sw.Old; % value is the current value of the variable including the boundary cells
% p.Old = cellBoundary(MeshStructure, BCp, p0);
p.Old=p0;
p.value = p.Old;
pgrad = gradientTerm(MeshStructure, p.value);
sw_ave = tvdMean(MeshStructure, sw.value,-pgrad, FL);
% pgrad = gradientTerm(MeshStructure, p.value);
Lw_ave = harmonicMean(MeshStructure, perm./mu_water);
Lo_ave = harmonicMean(MeshStructure, perm./mu_oil);
% explicit source term
s_ex = zeros(m+2,n+2);
s_ex(2:m+1,2:n+1) = s;
RHSs = s_ex(:); % explicit source term to be added to the rhs
% start the loop
error1 = 1e5;
while error1>eps1
    % step 1) calculate the average values
    sw_ave = tvdMean(MeshStructure, sw.value, -pgrad, FL);
    Lw = Lw_ave.*krw(sw_ave);
    Lo = Lo_ave.*kro(sw_ave);

    % step 2) calculate the pressure profile
    L_ave = Lo+Lw;
    Meq = diffusionTerm(MeshStructure, L_ave);
    [BCMp, BCRHSp] = boundaryCondition(MeshStructure, BCp);

    % solve the linear system of equations and reshape the result
    Mp = Meq + BCMp;
    RHSp = BCRHSp - RHSs; % the whole continuity is multiplied by a minus sign
    P = Mp\RHSp;
    p.value = reshape(full(P), m+2, n+2);

    pgrad = gradientTerm(MeshStructure, p.value);
    for j = 1:3
        sw_ave = tvdMean(MeshStructure, sw.value,-pgrad, FL);
    %     sw_ave = arithmeticMean(MeshStructure, sw.value);
        Lw = Lw_ave.*krw(sw_ave);
    %     Lo = Lo_ave.*kro(sw_ave);
        % step 3) calculate the new value of sw
        [Mtrans, RHStrans] = transientTerm(MeshStructure, ones(m,n), dt, sw);
        u = -dkrwdsw(sw_ave).*Lw_ave.*pgrad;
        [Mconv, RHSconv] = convectionTvdTerm(MeshStructure, u, sw.value, FL);
    %     Mconv = convectionTerm(MeshStructure, u);
        facevar = (-Lw+dkrwdsw(sw_ave).*Lw_ave.*sw_ave).*pgrad;
        RHSdiv = divergenceTerm(MeshStructure, facevar);
        [BCM, BCRHS] = boundaryCondition(MeshStructure, BCs);

        % construct the linear system
        M = Mtrans+Mconv+BCM;
        RHS = RHStrans+BCRHS-RHSdiv+RHSs+RHSconv;
        SW = M\RHS;
        sw_new = reshape(full(SW), m+2, n+2);
        error1 = sum(sum(abs(sw_new-sw.value)));
        sw.value = sw_new;
    end
end
end

% define the relperm curves


% function r = krw(sw)
% r.xvalue = sw.xvalue.^4;
% r.yvalue = sw.yvalue.^4;
% end
%
% function r = kro(sw)
% r.xvalue = (1-sw.xvalue.^2).*(1-sw.xvalue).^2;
% r.yvalue = (1-sw.yvalue.^2).*(1-sw.yvalue).^2;
% end
%
% function r = dkrwdsw(sw)
% r.xvalue = 4*sw.xvalue.^3;
% r.yvalue = 4*sw.yvalue.^3;
% end


function r = krw(sw)
fmmob = 10000;
fmdry = 0.2 ;
epdry = 1000 ;

fm.xvalue =  (sw.xvalue>fmdry).*(1 + fmmob.*atan(epdry.*(sw.xvalue-fmdry))./pi)+(sw.xvalue<=fmdry);
fm.yvalue =  (sw.yvalue>fmdry).*(1 + fmmob.*atan(epdry.*(sw.yvalue-fmdry))./pi)+(sw.yvalue<=fmdry);
r.xvalue = sw.xvalue.^2./fm.xvalue;
r.yvalue = sw.yvalue.^2./fm.yvalue;
end

function r = kro(sw)
r.xvalue = (1-sw.xvalue.^4).*(1-sw.xvalue).^4;
r.yvalue = (1-sw.yvalue.^4).*(1-sw.yvalue).^4;
end

function r = dkrwdsw(sw)
fmmob = 10000;
fmdry = 0.2 ;
epdry = 1000 ;
fm.xvalue =  (sw.xvalue>fmdry).*(1 + fmmob.*atan(epdry.*(fmdry-sw.xvalue))./pi)+(sw.xvalue<=fmdry);
fm.yvalue =  (sw.yvalue>fmdry).*(1 + fmmob.*atan(epdry.*(fmdry-sw.yvalue))./pi)+(sw.yvalue<=fmdry);
r.xvalue = 2*sw.xvalue./fm.xvalue;
r.yvalue = 2*sw.yvalue./fm.yvalue;
end
