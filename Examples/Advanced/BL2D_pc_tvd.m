% 2D Buckley-Leverett solution including source and  
% Geometry
m = 60; % number of cells in x direction
n = 20; % number of cells in y direction
W = 50; % width
H = 20; % height
MeshStructure = buildMesh2D(m, n, W, H);


% initial condition
so_init = ones(m,n);
p0 = 100e5;
p_init = p0*ones(m,n);

% rock parameters
k = 1e-12; % [m^2]
phi = 0.3; % [-] average porosity
C_DP = 0.7; % Dijkstra Parsons coef
[perm, poros] = RandPermField(k, phi, C_DP, 30, 10, m, n);
relperm_opts = relpermOptions();
pc_opts = capillaryPressureOptions();

% fluid properties
mu_water = ones(m,n)*1e-3; % water viscosity
mu_oil = ones(m,n)*30e-3; % oil viscosity
rho_oil = 800; % [kg/m^3]
rho_water = 1000; % [kg/m^3]

% boundary conditions and source terms
BCp = createBC(MeshStructure); % all Neumann BC for pressure
BCs = BCp; % saturation BC
n_bc = 3;
% BCp.left(1:n_bc) = 1; % change to Dirichlet
BCs.left.a(1:n_bc) = 0;  BCs.left.b(1:n_bc) = 1; BCs.left.c(1:n_bc) = 0;

qo = zeros(m,n);
qw = zeros(m,n);
% qw(1:2,1:2) = 1e-4;

% Define some of the variables over the domain
% gravity term:
g_val= 0;%9.81; % m/s^2
g.xvalue = zeros(m+1,n);
g.yvalue = -ones(m,n+1)*g_val; % m/s^2

% Mobility terms
Lw_face = harmonicMean(MeshStructure, perm./mu_water);
Lo_face = harmonicMean(MeshStructure, perm./mu_oil);

% define relperm functions
relperm_opts.phase = 'oil';
pc_opts.phase.phase = 'oil';
pc_opts.OilWaterExponent = 2;

% initialize
% so.Old = cellBoundary(MeshStructure, BCs, so_init);
so.Old = ones(m+2,n+2);
so.value = so.Old; % value is the current value of the variable including the boundary cells
so_k = so.value;
% p.Old = cellBoundary(MeshStructure, BCp, p_init);
p.Old = p0*ones(m+2,n+2);
p.value = p.Old;
u = -(gradientTerm(MeshStructure, p.value)-rho_oil*g);
% so_face = tvdMean(MeshStructure, so.value, u, FL);
so_face = arithmeticMean(MeshStructure, so.value);
% explicit source term
RHS_qo = sourceExplicitTerm(MeshStructure, qo);
RHS_qw = sourceExplicitTerm(MeshStructure, qw);
FL = fluxLimiter('MUSCL');
% start the loop
eps1 = 1e-5;
dt = 500;
t_end = 1000*dt;
for t = 0:dt:t_end
    error2 = 1e5;
    while error2>eps1
        % step 1) calculate the average values
%         so_face = arithmeticMean(MeshStructure, so.value);
        [krw, ~] = wateroilrelperm(so_face, relperm_opts);
        [kro, dkro] = oilwaterrelperm(so_face, relperm_opts);
        [~, dpc] = wateroilcapillarypressure(so_face, pc_opts);
        pcgrad = dpc.*gradientTerm(MeshStructure, so.value);

        Lw = Lw_face.*krw;
        Lo = Lo_face.*kro;
        L_rho_face = -(rho_oil*Lo+rho_water*Lw).*g+Lo.*pcgrad;
        L = Lo+Lw;
        % update the no flow boundary condition
        BCp.top.c(:) = -L_rho_face.yvalue(:,end)./L.yvalue(:,end);
        BCp.bottom.c(:) = -L_rho_face.yvalue(:,1)./L.yvalue(:,1);
        BCp.left.c(:) = -L_rho_face.xvalue(1,:)./L.xvalue(1,:);
        BCp.right.c(:) = -L_rho_face.xvalue(end,:)./L.xvalue(end,:);
        BCp.left.a(1:n_bc) = 0;  BCp.left.b(1:n_bc) = 1; BCp.left.c(1:n_bc) = 10*p0; % Pa, Dirichlet
        BCp.right.a(n-n_bc:n) = 0;  BCp.right.b(n-n_bc:n) = 1; BCp.right.c(n-n_bc:n) = p0; % Pa, Dirichlet

        % step 2) calculate the pressure profile
        Mp = diffusionTerm(MeshStructure, -L);
        [BCMp, BCRHSp] = boundaryCondition(MeshStructure, BCp);
        RHSp = divergenceTerm(MeshStructure, L_rho_face);
        % solve the linear system of equations and reshape the result
        Mpt = Mp + BCMp;
        RHSpt = BCRHSp + RHSp + RHS_qo + RHS_qw; % the whole continuity is multiplied by a minus sign
        P = Mpt\RHSpt;
        p.value = reshape(full(P), m+2, n+2);
        pgrad = (gradientTerm(MeshStructure, p.value)-rho_oil*g);
        error1 = 1e5;
        for j = 1:5
            % step 3) calculate the new value of sw
            so_face = tvdMean(MeshStructure, so.value, u, FL);
            [krw, ~] = wateroilrelperm(so_face, relperm_opts);
            [kro, dkro] = oilwaterrelperm(so_face, relperm_opts);
            [~, dpc] = wateroilcapillarypressure(so_face, pc_opts);
            pcgrad = dpc.*gradientTerm(MeshStructure, so.value);

            Lw = Lw_face.*krw;
            Lo = Lo_face.*kro;
            
            [Mtrans, RHStrans] = transientTerm(MeshStructure, poros, dt, so);
            u = -dkro.*Lo_face.*pgrad;
            [Mconv, RHSconv] = convectionTvdTerm(MeshStructure, u, so.value, FL);
            Mdiff = diffusionTerm(MeshStructure, -Lo.*dpc);
            facevar = (-Lo+dkro.*Lo_face.*so_face).*pgrad;
            RHSdiv = divergenceTerm(MeshStructure, facevar);
            [BCM, BCRHS] = boundaryCondition(MeshStructure, BCs);

            % construct the linear system
            M = Mtrans+Mconv+Mdiff+BCM;
            RHS = RHStrans+BCRHS-RHSdiv+RHS_qo+RHSconv;
            SO = M\RHS;
%             SO = agmg(M,RHS, [], 1e-7, 200);
            so_new = reshape(full(SO), m+2, n+2);
            error1 = sum(sum(abs(so_new-so.value)));
            so.value = so_new;
  
        end
        error2 = sum(sum(abs(so_new-so_k)));
        so_k = so_new;
        disp(error2);
    end
    so.Old = so.value;
    figure(1); pcolor(so.value(2:end-1,2:end-1)'); colorbar; shading flat
end