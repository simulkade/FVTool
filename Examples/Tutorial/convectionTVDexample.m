% another simple 1D convection example with periodic BC
% diffusion coefficients are defined but not used
% It compares the results of upwind with TVD and shows how diffusive the upwind
% scheme can be
clc; clear;
% define a 1D domain and mesh
W = 1;
Nx = 500;
mesh1 = createMesh1D(Nx, W);
x = mesh1.cellcenters.x;
% define the boundaries
BC = createBC(mesh1); % all Neumann
BC.left.periodic=1;
BC.right.periodic =1;
% Initial values
phi_old = createCellVariable(mesh1, 0.0, BC);
phi_old.value(20:120) = 1;
phi_old.value(180:400)= sin(x(180:400)*10*pi());
% initial guess for phi
phi = phi_old;
% initial values for upwind scheme
phiuw = phi;
phiuw_old=phi;
% keep the initial values for visualization
phiinit=phi_old;
% velocity field
u = 0.3;
uf = createFaceVariable(mesh1, u);
% diffusion field
D = 1e-2;
Df = createFaceVariable(mesh1, D);
% transient term coefficient
alfa = createCellVariable(mesh1,1.0);
% upwind convection term
Mconvuw = convectionUpwindTerm1D(uf);
% define the BC term
[Mbc, RHSbc] = boundaryCondition(BC);
% choose a flux limiter
FL = fluxLimiter('Superbee');
% solver
dt = 0.0005; % time step
final_t = W/u;
t = 0;
while t<final_t
    t = t+dt;
    % inner loop for TVD scheme
    [Mt, RHSt] = transientTerm(phi_old, dt, alfa);
    for j = 1:5
        [Mconv, RHSconv] = convectionTvdTerm1D(uf, phi, FL);
        M = Mconv+Mt+Mbc;
        RHS = RHSt+RHSbc+RHSconv;
        phi = solvePDE(mesh1, M, RHS);
    end
    [Mtuw, RHStuw] = transientTerm(phiuw_old, dt, alfa);
    Muw = Mconvuw+Mtuw+Mbc;
    RHSuw = RHStuw+RHSbc;
    phiuw = solvePDE(mesh1, Muw, RHSuw);
    phiuw_old = phiuw;
    phi_old = phi;
    figure(1);plot(x, phiinit.value(2:Nx+1), x, phi.value(2:Nx+1), '-o', x, ...
        phiuw.value(2:Nx+1)); drawnow;

end
