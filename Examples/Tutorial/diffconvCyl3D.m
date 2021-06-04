%% Transient diffusion equation
%% PDE and boundary conditions
% The transient diffusion equation reads
%
% $$\alpha\frac{\partial c}{\partial t}+\nabla.\left(-D\nabla c\right)=0,$$
%
% where $c$ is the independent variable (concentration, temperature, etc)
% , $D$ is the diffusion coefficient, and $\alpha$ is a constant.
% This one is cylindrical; you can plot it better if you want by modifying
% the visualization function
% Written by Ali A. Eftekhari
% Last checked: June 2021
clc

%% Define the domain and create a mesh structure
L = 50;  % domain length
Nx = 10; % number of cells
m = createMeshCylindrical3D(Nx,Nx+3,Nx+6, L/2,2*pi,L);
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
u = createFaceVariable(m, [0,0,0.1]);
%% define initial values
c_init = 1;
c_old = createCellVariable(m, c_init,BC); % initial values
c = c_old; % assign the old value of the cells to the current values
%% loop
Dave = harmonicMean(D);
Mdiff = diffusionTerm(Dave);
[Mbc, RHSbc] = boundaryCondition(BC);
FL = fluxLimiter('Superbee');
Mconv = convectionTvdTerm(u, c, FL);
dt = 1; % time step
final_t = 50;
for t=dt:dt:final_t
    [M_trans, RHS_trans] = transientTerm(c, dt, alfa);
    M = M_trans-Mdiff+Mbc+Mconv;
    RHS = RHS_trans+RHSbc;
    c = solvePDE(m,M, RHS);
    c_old = c;
end
%% visualization
 figure(1);visualizeCells(c);
