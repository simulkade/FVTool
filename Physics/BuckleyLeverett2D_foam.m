function [sg, p] = BuckleyLeverett2D_foam(MeshStructure, BCp, BCs, ...
			perm, mu_water, mu_gas, sg0, p0, dt, s, eps1)
% This function solves the continuity equation for an incompressible fluid 
% flow in a porous medium, where Darcy's law is applicable.
 
% extract the mesh data
mn = MeshStructure.numberofcells;
m = mn(1); n = mn(2);



% start the solver here
% first solve the convection equation for sw
% initialize
% sw.Old = cellBoundary(MeshStructure, BCs, sw0);
sg.Old = sg0;
sg.value = sg.Old; % value is the current value of the variable including the boundary cells
% p.Old = cellBoundary(MeshStructure, BCp, p0);
p.Old=p0;
p.value = p.Old;
% pgrad = gradientTerm(MeshStructure, p.value);
Lg_ave = harmonicMean(MeshStructure, perm./mu_gas);
Lw_ave = harmonicMean(MeshStructure, perm./mu_water);
% explicit source term
s_ex = zeros(m+2,n+2);
s_ex(2:m+1,2:n+1) = s;
RHSs = s_ex(:); % explicit source term to be added to the rhs
% start the loop
error1 = 1e5;
while error1>eps1
    % step 1) calculate the average values
    sg_ave = arithmeticMean(MeshStructure, sg.value);
    Lg = Lg_ave.*krg(sg_ave);
    Lw = Lw_ave.*krw(sg_ave);
    
    % step 2) calculate the pressure profile
    L_ave = Lw+Lg;
    Meq = diffusionTerm(MeshStructure, L_ave);
    [BCMp, BCRHSp] = boundaryCondition(MeshStructure, BCp);

    % solve the linear system of equations and reshape the result
    Mp = Meq + BCMp;
    RHSp = BCRHSp - RHSs; % the whole continuity is multiplied by a minus sign
    P = Mp\RHSp;
    p.value = reshape(full(P), m+2, n+2);
    
    pgrad = gradientTerm(MeshStructure, p.value);
    sg_ave = upwindMean(MeshStructure, -pgrad, sg.value);
%     sw_ave = arithmeticMean(MeshStructure, sw.value);
    Lg = Lg_ave.*krg(sg_ave);
%     Lo = Lo_ave.*kro(sw_ave);
    % step 3) calculate the new value of sw
    [Mtrans, RHStrans] = transientTerm(MeshStructure, ones(m,n), dt, sg);
    u = -dkrgdsg(sg_ave).*Lg_ave.*pgrad;
    Mconv = convectionUpwindTerm(MeshStructure, u);
%     Mconv = convectionTerm(MeshStructure, u);
    facevar = (-Lg+dkrgdsg(sg_ave).*Lg_ave.*sg_ave).*pgrad;
    RHSdiv = divergenceTerm(MeshStructure, facevar);
    [BCM, BCRHS] = boundaryCondition(MeshStructure, BCs);

    % construct the linear system
    M = Mtrans+Mconv+BCM;
    RHS = RHStrans+BCRHS-RHSdiv+RHSs;
    SG = M\RHS;
    sg_new = reshape(full(SG), m+2, n+2);
    error1 = sum(sum(abs(sg_new-sg.value)))
    sg.value = sg_new;
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


function r = krg(sg)
fmmob = 100000;
fmdry = 0.2 ;
epdry = 10 ;
swc = 0;
sgr = 0;
krg0 = 1;
ng = 2;
kr = krg0*(1-(1-sg-swc)/(1-sgr-swc)).^ng;
fm = (1-sg>fmdry).*(1+fmmob*(0.5+atan(epdry.*(1-sg-fmdry))/pi()))+(1-sg<=fmdry);
r= kr./fm;
end

function r = krw(sg)
krw0 = 1;
nw = 4;
sgr = 0;
swc = 0;
r = krw0*((1-sg-swc)/(1-sgr-swc)).^nw;
end

function r = dkrgdsg(sg)
fmmob = 100000;
fmdry = 0.2 ;
epdry = 10 ;
swc = 0;
sgr = 0;
krg0 = 1;
ng = 2;
kr = krg0*(1-(1-sg-swc)/(1-sgr-swc)).^ng;
dkrdsg = (krg0*ng*(1-(-swc-sg+1)/(-swc-sgr+1)).^(ng-1))/(-swc-sgr+1);
fm = (1-sg>fmdry).*(1+fmmob*(0.5+atan(epdry.*(1-sg-fmdry))/pi()))+(1-sg<=fmdry);
dfmdsg = (1-sg>fmdry).*(-(epdry*fmmob)./(pi*(epdry^2*(-sg-fmdry+1).^2+1)));
r = (dkrdsg.*fm-dfmdsg.*kr)./fm.^2;
end
