% a tutorial adapted from the fipy convection diffusion 1D example
% see: http://www.ctcms.nist.gov/fipy/examples/convection/index.html

clc
clear
%% define the domain
L = 1;  % domain length
Nx = 25; % number of cells
meshstruct = buildMesh1D(Nx, L);
BC = createBC(meshstruct); % all Neumann boundary condition structure
BC.left.a = 0; BC.left.b=1; BC.left.c=0; % left boundary
BC.right.a = 0; BC.right.b=1; BC.right.c=1; % right boundary
x = meshstruct.cellcenters.x;
%% define the transfer coeffs
D_val = -1;
D = createCellVariable(meshstruct, D_val);
Dave = harmonicMean(meshstruct, D); % convert a cell variable to face variable
alfa = createCellVariable(meshstruct, 1);
u = -10;
u_face = createFaceVariable(meshstruct, u);
%% solve
Mconv =  convectionTerm(meshstruct, u_face);
Mconvupwind =  convectionUpwindTerm(meshstruct, u_face);
Mdiff = diffusionTerm(meshstruct, Dave);
[Mbc, RHSbc] = boundaryCondition(meshstruct, BC);
M = Mconv-Mdiff-Mbc;
Mupwind = Mconvupwind-Mdiff-Mbc;
RHS = -RHSbc;
c = solvePDE(meshstruct, M, RHS);
c_upwind = solvePDE(meshstruct, Mupwind, RHS);
c_analytical = (1-exp(u*x/D_val))/(1-exp(u*L/D_val));
figure(1);plot(x, c(2:Nx+1), x, c_upwind(2:Nx+1), '-r', x, c_analytical, '^');
legend('central', 'upwind', 'analytical');



