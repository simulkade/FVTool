function [sw, p] = BuckleyLeverett1D(MeshStructure, BCp, BCs, ...
			perm, poros, mu_oil, mu_water, sw0, p0, dt, eps1)
% This function solves the continuity equation for an incompressible fluid
% flow in a porous medium, where Darcy's law is applicable.
% Copyright (c) 2012-2016 Ali Akbar Eftekhari
% See the license file

% extract the mesh data
m = MeshStructure.numberofcells;

% define the relperm curves
krw = @(sw)(sw.^4);
kro = @(sw)((1-sw.^2).*(1-sw).^2);
dkrwdsw = @(sw)(4*sw.^3);

% start the solver here
% first solve the convection equation for sw
% initialize
% sw.Old = cellBoundary(MeshStructure, BCs, sw0);
sw.Old = sw0;
sw.value = sw.Old; % value is the current value of the variable including the boundary cells
% p.Old = cellBoundary(MeshStructure, BCp, p0);
p.Old = p0;
p.value = p.Old;
pgrad = gradientTerm(MeshStructure, p.value);
% start the loop
error1 = 1e5;
while error1>eps1
    % step 1) calculate the average values
    sw_ave = upwindMean(MeshStructure, -pgrad, sw.value);
    Lw_ave = harmonicMean(MeshStructure, perm./mu_water);
    Lw = Lw_ave.xvalue.*krw(sw_ave.xvalue);
    Lo_ave = harmonicMean(MeshStructure, perm./mu_oil);
    Lo = Lo_ave.xvalue.*kro(sw_ave.xvalue);

    % step 2) calculate the pressure profile
    L_ave.xvalue = Lo+Lw;
    Meq = diffusionTerm(MeshStructure, L_ave);
    [BCMp, BCRHSp] = boundaryCondition(MeshStructure, BCp);

    % solve the linear system of equations and reshape the result
    Mp = Meq + BCMp;
    RHSp = BCRHSp;
    P = Mp\RHSp;
    p.value = reshape(full(P), m+2, 1);

    pgrad = gradientTerm(MeshStructure, p.value);
    sw_ave = upwindMean(MeshStructure, -pgrad, sw.value);
    Lw_ave = harmonicMean(MeshStructure, perm./mu_water);
    Lw = Lw_ave.xvalue.*krw(sw_ave.xvalue);
%     u.xvalue = -(Lo+Lw).*pgrad.xvalue;

    % step 3) calculate the new value of sw
    [Mtrans, RHStrans] = transientTerm(MeshStructure, poros, dt, sw);
    u_temp.xvalue = -dkrwdsw(sw_ave.xvalue).*Lw_ave.xvalue.*pgrad.xvalue;
    Mconv = convectionUpwindTerm(MeshStructure, u_temp);
    facevar.xvalue = (-Lw+dkrwdsw(sw_ave.xvalue).*Lw_ave.xvalue.*sw_ave.xvalue).*pgrad.xvalue;
    RHSdiv = divergenceTerm(MeshStructure, facevar);
    [BCM, BCRHS] = boundaryCondition(MeshStructure, BCs);

    % construct the linear system
    M = Mtrans+Mconv+BCM;
    RHS = RHStrans+BCRHS-RHSdiv;
    SW = M\RHS;
    sw_new = reshape(full(SW), m+2, 1);
    error1 = sum(abs(sw_new-sw.value));
    sw.value = sw_new;
end
