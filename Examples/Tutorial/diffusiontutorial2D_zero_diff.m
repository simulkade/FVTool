%% Transient diffusion equation
%% PDE and boundary conditions
% The transient diffusion equation reads
%
% $$\alpha\frac{\partial c}{\partial t}+\nabla.\left(-D\nabla c\right)=0,$$
%
% where $c$ is the independent variable (concentration, temperature, etc)
% , $D$ is the diffusion coefficient, and $\alpha$ is a constant.
% Written by Ali A. Eftekhari
% Last checked: June 2021
clc

%% Define the domain and create a mesh structure
L = 50;  % domain length
Nx = 20; % number of cells
n_zero = 5:8;
m = createMesh2D(Nx,Nx, L,L);
%% Create the boundary condition structure
BC = createBC(m); % all Neumann boundary condition structure
BC.left.a(:) = 0; BC.left.b(:)=1; BC.left.c(:)=0; % left boundary
BC.right.a(:) = 0; BC.right.b(:)=1; BC.right.c(:)=0; % right boundary
BC.top.a(:) = 0; BC.top.b(:)=1; BC.top.c(:)=0; % top boundary
BC.bottom.a(:) = 0; BC.bottom.b(:)=1; BC.bottom.c(:)=0; % bottom boundary
%% define the transfer coeffs
D_val = 1;
D = createCellVariable(m, D_val);
D.value(n_zero, n_zero) = 0.0; % zero diffusion; no flux into/from cells
Dave = harmonicMean(D);
Dave.xvalue(isnan(Dave.xvalue)) = 0.0;
Dave.yvalue(isnan(Dave.yvalue)) = 0.0;
alfa = createCellVariable(m, 1);
%% define initial values
c_init = 1;
c_old = createCellVariable(m, c_init,BC); % initial values
c = c_old; % assign the old value of the cells to the current values
%% loop
dt = 1; % time step
final_t = 100;
for t=dt:dt:final_t
    [M_trans, RHS_trans] = transientTerm(c_old, dt, alfa);
    Mdiff = diffusionTerm(Dave);
    [Mbc, RHSbc] = boundaryCondition(BC);
    M = M_trans-Mdiff+Mbc;
    RHS = RHS_trans+RHSbc;
    c = solvePDE(m,M, RHS);
    c_old = c;
    figure(1);visualizeCells(c);drawnow;
end
