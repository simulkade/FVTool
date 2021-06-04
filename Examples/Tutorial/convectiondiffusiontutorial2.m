% a tutorial adapted from the fipy convection diffusion 1D example
% see: http://www.ctcms.nist.gov/fipy/examples/convection/index.html
% Written by Ali A. Eftekhari
% Last checked: June 2021
clc
%% define the domain
L = 1;  % domain length
Nx = 25; % number of cells
meshstruct = createMesh1D(Nx, L);
BC = createBC(meshstruct); % all Neumann boundary condition structure
BC.left.a = 0; BC.left.b=1; BC.left.c=1; % left boundary
% BC.right.a = 0; BC.right.b=1; BC.right.c=1; % right boundary
x = meshstruct.cellcenters.x;
%% define the transfer coeffs
D_val = 1e-8;
D = createCellVariable(meshstruct, D_val);
Dave = harmonicMean(D); % convert a cell variable to face variable
u = 7.6e-6;
u_face = createFaceVariable(meshstruct, u);
alfa = createCellVariable(meshstruct,1.0);
phi_old = createCellVariable(meshstruct, 0.0, BC);
phi_old_upwind = phi_old;
phi_tvd = phi_old;
phi_old_tvd = phi_old;
%% solve
dt = 10; % time step
final_t = 100000;
t = 0;
FL = fluxLimiter('Superbee');

Mconv =  convectionTerm(u_face);
Mconvupwind =  convectionUpwindTerm(u_face);
Mdiff = diffusionTerm(Dave);
[Mbc, RHSbc] = boundaryCondition(BC);
    
while t<final_t
    disp(t)
    t=t+dt;
    [Mt, RHSt] = transientTerm(phi_old, dt, alfa);
    [Mt_tvd, RHSt_tvd] = transientTerm(phi_old_tvd, dt, alfa);
    [Mt_upwind, RHSt_upwind] = transientTerm(phi_old_upwind, dt, alfa);
    M = Mconv-Mdiff+Mbc+Mt;
    Mupwind = Mconvupwind-Mdiff+Mbc+Mt_upwind;
    RHS = RHSbc+RHSt;
    RHS_upwind = RHSbc+RHSt_upwind;
    phi = solvePDE(meshstruct, M, RHS);
    phi_upwind = solvePDE(meshstruct, Mupwind, RHS_upwind);
    for j = 1:5
        [Mconv, RHSconv] = convectionTvdTerm1D(u_face, phi_tvd, FL);
        M_tvd = Mconv+Mt_tvd+Mbc;
        RHS_tvd = RHSt_tvd+RHSbc+RHSconv;
        phi_tvd = solvePDE(meshstruct, M_tvd, RHS_tvd);
    end
    phi_old = phi;
    phi_old_upwind = phi_upwind;
    phi_old_tvd = phi_tvd;
end
visualizeCells(phi);
hold on
visualizeCells(phi_upwind);
hold on
visualizeCells(phi_tvd);
hold off
% 
%     figure(1);plot(x, phi.value(2:Nx+1), 'linewidth', 2, ...
%      x, phi_upwind.value(2:Nx+1), '-r', 'linewidth', 2)
%     legend('central', 'upwind', 'analytical');
