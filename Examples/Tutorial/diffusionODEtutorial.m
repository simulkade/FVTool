% a tutorial adapted from the fipy diffusion 1D example
% see: http://www.ctcms.nist.gov/fipy/examples/diffusion/index.html
% How to convert a PDE into an ODE and solve it using a matlab ode solver
% needs more documentation
clc
clear

%% define the domain
L = 50;  % domain length
Nx = 20; % number of cells
meshstruct = createMesh1D(Nx, L);
BC = createBC(meshstruct); % all Neumann boundary condition structure
BC.left.a = 0; BC.left.b=1; BC.left.c=1; % left boundary
BC.right.a = 0; BC.right.b=1; BC.right.c=0; % right boundary
x = meshstruct.cellcenters.x;
%% define the transfer coeffs
D_val = 1;
D = createFaceVariable(meshstruct, D_val);
alfa = createCellVariable(meshstruct, 1);
%% define initial values
c_old = createCellVariable(meshstruct, 0); % initial values
c_value = c_old;
%% loop
dt = 0.1; % time step
final_t = 100;
Mdiff = diffusionTerm(D);
[M, RHS] = combineBC1D(meshstruct, BC, Mdiff, ...
    zeros(Nx+2,1));
% eq: dcdt = D d2c/dx2
dcdt = @(t,c)(M*c-RHS);
[t_temp, c_temp] = ode45(dcdt, [0 final_t], internalCells(c_old));
c_analytical = 1-erf(x./(2*sqrt(D_val*final_t)));
plot(x, c_temp(end,:), x, c_analytical, 'o')