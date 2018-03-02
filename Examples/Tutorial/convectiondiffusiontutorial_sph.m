% convection diffusion in a spherical coodinate (1D)

clc
clear
%% define the domain
L = 1;  % domain length
Nx = 25; % number of cells
meshstruct = createMeshSpherical1D(Nx, L);
meshstruct.cellcenters.x=meshstruct.cellcenters.x+0.1;
meshstruct.facecenters.x=meshstruct.facecenters.x+0.1;

BC = createBC(meshstruct); % all Neumann boundary condition structure
BC.left.a = 0; BC.left.b=1; BC.left.c=0; % left boundary
BC.right.a = 0; BC.right.b=1; BC.right.c=1; % right boundary
x = meshstruct.cellcenters.x;

%% define the transfer coeffs
D_val = 1;
D = createCellVariable(meshstruct, D_val);
Dave = harmonicMean(D); % convert a cell variable to face variable
alfa = createCellVariable(meshstruct, 1);
u = -5;
u_face = createFaceVariable(meshstruct, u);
%% solve
Mconv =  convectionTerm(u_face);
Mconvupwind =  convectionUpwindTerm(u_face);
Mdiff = diffusionTerm(Dave);
[Mbc, RHSbc] = boundaryCondition(BC);
M = Mconv-Mdiff-Mbc;
Mupwind = Mconvupwind-Mdiff-Mbc;
RHS = -RHSbc;
c = solvePDE(meshstruct, M, RHS);
c_upwind = solvePDE(meshstruct, Mupwind, RHS);
% c_analytical = (1-exp(u*x/D_val))/(1-exp(u*L/D_val));
figure(1);
plot(x, c.value(2:Nx+1), ... % 'linewidth', 2, ...
 x, c_upwind.value(2:Nx+1), '-r', 'linewidth', 2)
%   x, c_analytical, '^', 'linewidth', 2);
legend('central', 'upwind');
