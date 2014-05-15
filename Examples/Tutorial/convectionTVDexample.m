% another simple 1D convection example with periodic BC
% diffusion coefficients are defined but not used
% It compares the results of upwind with TVD and shows how diffusive the upwind
% scheme can be
clc; clear;
% define a 1D domain and mesh
W = 1;
Nx = 500;
mesh1 = buildMesh1D(Nx, W);
x = mesh1.cellcenters.x;
% define the boundaries
BC = createBC(mesh1); % all Neumann
BC.left.periodic=1;
BC.right.periodic =1;
% Initial values
phi.Old = createCellVariable(mesh1, 0.0, BC);
phi.Old(20:120) = 1;
phi.Old(180:400)= sin(x(180:400)*10*pi());
% initial guess for phi
phi.value = phi.Old;
% initial values for upwind scheme
phiuw = phi;
% keep the initial values for visualization
phiinit=phi.Old;
% velocity field
u = 0.3;
uf = createFaceVariable(mesh1, u);
% diffusion field
D = 1e-2;
Df = createFaceVariable(mesh1, D);
% transient term coefficient
alfa = createCellVariable(mesh1,1.0);
% upwind convection term
Mconvuw = convectionUpwindTerm1D(mesh1, uf);
% define the BC term
[Mbc, RHSbc] = boundaryCondition(mesh1, BC);
% choose a flux limiter
FL = fluxLimiter('MUSCL');
% solver
dt = 0.0005; % time step
final_t = W/u;
t = 0;
while t<final_t
    t = t+dt;
    % inner loop for TVD scheme
    for j = 1:5
        [Mt, RHSt] = transientTerm(mesh1, alfa, dt, phi);
        [Mconv, RHSconv] = convectionTvdTerm1D(mesh1, uf, phi.value, FL);
        M = Mconv+Mt+Mbc;
        RHS = RHSt+RHSbc+RHSconv;
        PHI = M\RHS;
        phi.value = reshape(PHI, Nx+2, 1);
    end
    [Mtuw, RHStuw] = transientTerm(mesh1, alfa, dt, phiuw);
    Muw = Mconvuw+Mtuw+Mbc;
    RHSuw = RHStuw+RHSbc;
    PHI = Muw\RHSuw;
    phiuw.value = reshape(PHI, Nx+2, 1);
    phiuw.Old = phiuw.value;
    phi.Old = phi.value;
    figure(1);plot(x, phiinit(2:Nx+1), x, phi.value(2:Nx+1), '-o', x, ...
        phiuw.value(2:Nx+1)); drawnow;
    
end
