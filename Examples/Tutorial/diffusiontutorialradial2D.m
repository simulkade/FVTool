%% Transient diffusion equation
%% PDE and boundary conditions
% The transient diffusion equation reads
%
% $$\alpha\frac{\partial c}{\partial t}+\nabla.\left(-D\nabla c\right)=0,$$
%
% where $c$ is the independent variable (concentration, temperature, etc)
% , $D$ is the diffusion coefficient, and $\alpha$ is a constant. 
clc; clear;

%% Define the domain and create a mesh structure
L = 10;  % domain length
Nx = 10; % number of cells
Ntetta = 20;
m = buildMeshRadial2D(Nx,Ntetta, L,1*pi);
%% Create the boundary condition structure
BC = createBC(m); % all Neumann boundary condition structure
m.cellcenters.x=m.cellcenters.x+1;
m.facecenters.x=m.facecenters.x+1;
BC.left.a(:) = 0; BC.left.b(:)=1; BC.left.c(:)=2; % left boundary
BC.right.a(:) = 0; BC.right.b(:)=1; BC.right.c(:)=0; % right boundary
BC.right.a(1:floor(Ntetta/2)) = 0; BC.right.b(1:floor(Ntetta/2))=1; BC.right.c(1:floor(Ntetta/2))=3; % right boundary
% BC.top.periodic=1;
BC.top.a(:) = 0; BC.top.b(:)=1; BC.top.c(:)=0; % top boundary
BC.bottom.a(:) = 0; BC.bottom.b(:)=1; BC.bottom.c(:)=0; % bottom boundary
%% define the transfer coeffs
D_val = 1;
D = createCellVariable(m, D_val);
alfa = createCellVariable(m, 1);
%% define initial values
c_init = 0.1;
c.Old = createCellVariable(m, c_init,BC); % initial values
c.value = c.Old; % assign the old value of the cells to the current values
%% loop
dt = 0.1; % time step
final_t = 20;
Dave = harmonicMean(m, D);
Mdiff = diffusionTerm(m, Dave);
[Mbc, RHSbc] = boundaryCondition(m, BC);
for t=dt:dt:final_t
    [M_trans, RHS_trans] = transientTerm(m, alfa, dt, c);
    M = M_trans-Mdiff+Mbc;
    RHS = RHS_trans+RHSbc;
    c.value = solvePDE(m,M, RHS);
    c.Old = c.value;
    figure(1);visualizeCells(m, c.value);shading interp; drawnow;
end
