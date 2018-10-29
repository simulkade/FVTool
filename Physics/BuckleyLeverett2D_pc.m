function [so, p] = BuckleyLeverett2D_pc(MeshStructure, BCp, BCs, ...
			perm, poros, mu_oil, mu_water, rho_oil, rho_water, ...
            so_init, p_init, dt, qo, qw, relperm_opts, pc_opts, eps1)
% This function solves the continuity equation for an incompressible fluid
% flow in a porous medium, where Darcy's law is applicable.
% Copyright (c) 2012-2016 Ali Akbar Eftekhari
% See the license file 
% extract the mesh data
mn = MeshStructure.numberofcells;
m = mn(1); n = mn(2);

% Define some of the variables over the domain
% gravity term:
g_val= 9.81; % m/s^2
g.xvalue = zeros(m+1,n);
g.yvalue = -ones(m,n+1)*g_val; % m/s^2

% Mobility terms
Lw_face = harmonicMean(MeshStructure, perm./mu_water);
Lo_face = harmonicMean(MeshStructure, perm./mu_oil);

% define relperm functions
relperm_opts.phase = 'oil';
pc_opts.phase.phase = 'oil';

% initialize
so.Old = cellBoundary(MeshStructure, BCs, so_init);
so.value = sw.Old; % value is the current value of the variable including the boundary cells
p.Old = cellBoundary(MeshStructure, BCp, p_init);
p.value = p.Old;
pgrad = (gradientTerm(MeshStructure, p.value)-rho_oil*g);

% explicit source term
RHS_qo = sourceExplicitTerm(MeshStructure, qo);
RHS_qw = sourceExplicitTerm(MeshStructure, qw);

% start the loop
error1 = 1e5;
while error1>eps1
    % step 1) calculate the average values
    so_face = upwindMean(MeshStructure, -pgrad, so.value);
    [krw, ~] = wateroilrelperm(so_face, relperm_opts);
    [kro, dkro] = oilwaterrelperm(so_face, relperm_opts);
    [~, dpc] = wateroilcapillarypressure(so_face, pc_opts);
    pcgrad = dpc.*gradientTerm(MeshStructure, so.value);

    Lw = Lw_face.*krw;
    Lo = Lo_face.*kro;

    % step 2) calculate the pressure profile
    L = Lo+Lw;
    Mp = diffusionTerm(MeshStructure, -L);
    [BCMp, BCRHSp] = boundaryCondition(MeshStructure, BCp);
    RHSp = divergenceTerm(MeshStructure, ...
        -(rho_oil*Lo+rho_water*Lw).*g+Lo.*pcgrad);
    % solve the linear system of equations and reshape the result
    Mpt = Mp + BCMp;
    RHSpt = BCRHSp + RHSp + RHS_qo + RHS_qw; % the whole continuity is multiplied by a minus sign
    P = Mpt\RHSpt;
    p.value = reshape(full(P), m+2, n+2);

    pgrad = (gradientTerm(MeshStructure, p.value)-rho_oil*g);

    % step 3) calculate the new value of sw
    [Mtrans, RHStrans] = transientTerm(MeshStructure, poros, dt, so);
    u = -dkro.*Lo_face.*pgrad;
    Mconv = convectionUpwindTerm(MeshStructure, u);
    Mdiff = diffusionTerm(MeshStructure, -Lo.*dpc);
    facevar = (-Lo+dkro.*Lo_face.*so_face).*pgrad;
    RHSdiv = divergenceTerm(MeshStructure, facevar);
    [BCM, BCRHS] = boundaryCondition(MeshStructure, BCs);

    % construct the linear system
    M = Mtrans+Mconv+Mdiff+BCM;
    RHS = RHStrans+BCRHS-RHSdiv+RHS_qo;
    SO = M\RHS;
    so_new = reshape(full(SO), m+2, n+2);
    error1 = sum(sum(abs(so_new-so.value)));
    so.value = so_new;
end
end
