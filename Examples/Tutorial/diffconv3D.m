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
L = 50;  % domain length
Nx = 10; % number of cells
m = buildMesh3D(Nx,Nx+3,Nx+6, L/2,2*pi,L);
%% Create the boundary condition structure
BC = createBC(m); % all Neumann boundary condition structure
% BC.left.a(:) = 0; BC.left.b(:)=1; BC.left.c(:)=0; % left boundary
BC.right.a(:) = 0; BC.right.b(:)=1; BC.right.c(:)=0; % right boundary
% BC.top.a(:) = 0; BC.top.b(:)=1; BC.top.c(:)=3; % top boundary
% BC.bottom.a(:) = 0; BC.bottom.b(:)=1; BC.bottom.c(:)=0; % bottom boundary
BC.top.periodic=1;
BC.back.a(:) = 0; BC.back.b(:)=1; BC.back.c(:)=3; % back boundary
BC.front.a(:) = 0; BC.front.b(:)=1; BC.front.c(:)=0; % front boundary
%% define the transfer coeffs
D_val = 1;
D = createCellVariable(m, D_val);
alfa = createCellVariable(m, 1);
u = createFaceVariable(m, [0,0,0.5]);
%% define initial values
c_init = 1;
c.Old = createCellVariable(m, c_init,BC); % initial values
c.value = c.Old; % assign the old value of the cells to the current values
%% loop
Dave = harmonicMean(m, D);
Mdiff = diffusionTerm(m, Dave);
[Mbc, RHSbc] = boundaryCondition(m, BC);
FL = fluxLimiter('MUSCL');
Mconv = convectionTvdTerm(m, u, c.value, FL);
dt = 1; % time step
final_t = 50;
for t=dt:dt:final_t
    [M_trans, RHS_trans] = transientTerm(m, alfa, dt, c);
    M = M_trans-Mdiff+Mbc+Mconv;
    RHS = RHS_trans+RHSbc;
    c.value = solvePDE(m,M, RHS);
    c.Old = c.value;
    figure(1);visualizeCells(m, c.value);drawnow;
end
%% visualization
 figure(1);visualizeCells(m, c.value);