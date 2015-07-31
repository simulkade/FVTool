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
Nx = 20; % number of cells
Ntetta = 21;
m = createMeshRadial2D(Nx,Ntetta, L,2*pi);
%% Create the boundary condition structure
BC = createBC(m); % all Neumann boundary condition structure
m.cellcenters.x=m.cellcenters.x+1;
m.facecenters.x=m.facecenters.x+1;
BC.left.a(1:3:20) = 0; BC.left.b(1:3:20)=1; BC.left.c(1:3:20)=10; % left boundary
BC.right.a(:) = 0; BC.right.b(:)=1; BC.right.c(:)=0; % right boundary
% BC.right.a(1:floor(Ntetta/2)) = 0; BC.right.b(1:floor(Ntetta/2))=1; BC.right.c(1:floor(Ntetta/2))=3; % right boundary
BC.top.periodic=1;
% BC.top.a(:) = 0; BC.top.b(:)=1; BC.top.c(:)=0; % top boundary
% BC.bottom.a(:) = 0; BC.bottom.b(:)=1; BC.bottom.c(:)=0; % bottom boundary
%% define the transfer coeffs
D_val = 0.01;
D = createCellVariable(m, D_val);
alfa = createCellVariable(m, 1);
u_val=1;
u = createFaceVariable(m, [u_val, 0]);
%% define initial values
c_init = 0.1;
c_old = createCellVariable(m, c_init,BC); % initial values
c = c_old; % assign the old value of the cells to the current values
%% loop
dt = 0.1; % time step
final_t = 7;
Dave = harmonicMean(D);
Mdiff = diffusionTerm(Dave);
[Mbc, RHSbc] = boundaryCondition(BC);
FL = fluxLimiter('Superbee');
Mconv = convectionUpwindTerm(u);
for t=dt:dt:final_t
    [M_trans, RHS_trans] = transientTerm(c, dt, alfa);
    M = M_trans+Mconv-Mdiff+Mbc;
    RHS = RHS_trans+RHSbc;
    c = solvePDE(m,M, RHS);
    c_old = c.value;
    figure(1);visualizeCells(c);shading interp; drawnow;
end
