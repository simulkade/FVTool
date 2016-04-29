function [sw, p] = BuckleyLeverett3D(MeshStructure, BCp, BCs, ...
			perm, mu_oil, mu_water, sw0, p0, dt, s, eps1)
% This function solves the continuity equation for an incompressible fluid
% flow in a porous medium, where Darcy's law is applicable.
% Copyright (c) 2012-2016 Ali Akbar Eftekhari
% See the license file 
% extract the mesh data
mn = MeshStructure.numberofcells;
m = mn(1); n = mn(2); w = mn(3);



% start the solver here
% first solve the convection equation for sw
% initialize
% sw.Old = cellBoundary(MeshStructure, BCs, sw0);
sw.Old = sw0;
sw.value = sw.Old; % value is the current value of the variable including the boundary cells
% p.Old = cellBoundary(MeshStructure, BCp, p0);
p.Old=p0;
p.value = p.Old;
% pgrad = gradientTerm(MeshStructure, p.value);
Lw_ave = harmonicMean(MeshStructure, perm./mu_water);
Lo_ave = harmonicMean(MeshStructure, perm./mu_oil);
% explicit source term
s_ex = zeros(m+2,n+2, w+2);
s_ex(2:m+1,2:n+1, 2:w+1) = s;
RHSs = s_ex(:); % explicit source term to be added to the rhs
% start the loop
error1 = 1e5;
while error1>eps1
    % step 1) calculate the average values
    sw_ave = arithmeticMean(MeshStructure, sw.value);
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
    p.value = reshape(full(P), m+2, n+2, w+2);

    pgrad = gradientTerm(MeshStructure, p.value);
    sw_ave = upwindMean(MeshStructure, -pgrad, sw.value);
    Lw = Lw_ave.*krw(sw_ave);
%     Lo = Lo_ave.*kro(sw_ave);
    % step 3) calculate the new value of sw
    [Mtrans, RHStrans] = transientTerm(MeshStructure, ones(m,n,w), dt, sw);
    u = -dkrwdsw(sw_ave).*Lw_ave.*pgrad;
    Mconv = convectionUpwindTerm(MeshStructure, u);
    facevar = (-Lw+dkrwdsw(sw_ave).*Lw_ave.*sw_ave).*pgrad;
    RHSdiv = divergenceTerm(MeshStructure, facevar);
    [BCM, BCRHS] = boundaryCondition(MeshStructure, BCs);

    % construct the linear system
    M = Mtrans+Mconv+BCM;
    RHS = RHStrans+BCRHS-RHSdiv+RHSs;
    SW = M\RHS;
%     SW = agmg(M, RHS, 1e-5, 200);
    sw_new = reshape(full(SW), m+2, n+2, w+2);
    error1 = sum(sum(sum(abs(sw_new-sw.value))));
    sw.value = sw_new;
end
end

% define the relperm curves
function r = krw(sw)
r.xvalue = sw.xvalue.^4;
r.yvalue = sw.yvalue.^4;
r.zvalue = sw.zvalue.^4;
end

function r = kro(sw)
r.xvalue = (1-sw.xvalue.^2).*(1-sw.xvalue).^2;
r.yvalue = (1-sw.yvalue.^2).*(1-sw.yvalue).^2;
r.zvalue = (1-sw.zvalue.^2).*(1-sw.zvalue).^2;
end

function r = dkrwdsw(sw)
r.xvalue = 4*sw.xvalue.^3;
r.yvalue = 4*sw.yvalue.^3;
r.zvalue = 4*sw.zvalue.^3;
end
