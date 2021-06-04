%% Transient diffusion equation
%% PDE and boundary conditions
% The transient diffusion equation reads
%
% $$\alpha\frac{\partial c}{\partial t}+\nabla.\left(-D\nabla c\right)=f,$$
%
% where $c$ is the independent variable (concentration, temperature, etc)
% , $D$ is the diffusion coefficient, and $\alpha$ is a constant.
% OK this is a strange example; I vary the velocity field and the diffusion 
% coefficient in each cell, so it shows how you can do the same.
% Written by Ali A. Eftekhari
% Last checked: June 2021
clc
H=@(x)((x>0)+0.5*ones(size(x)).*(x==0));
%% Define the domain and create a mesh structure
L = 10;  % domain length
Nx = 20; % number of cells
m = createMesh1D(Nx,L);
X=m.cellcenters.x;
Xf=m.facecenters.x;
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
alfa = createCellVariable(m, 1);
u = createFaceVariable(m, 0.0);
%% define initial values
c_init = 1;
c_old = createCellVariable(m, c_init,BC); % initial values
c = c_old; % assign the old value of the cells to the current values
% assign different values of velocity and diffusivity
D = createCellVariable(m, 0.1+0.01*sin(X));
Dave = harmonicMean(D);
u.xvalue=0.5+0.03*cos(Xf); % they all shoud be on the cell faces
Mdiff = diffusionTerm(Dave);

f=createCellVariable(m,0.0); % source term
RHS_f=constantSourceTerm(f);
%% loop
[Mbc, RHSbc] = boundaryCondition(BC);
FL = fluxLimiter('Superbee');
dt = 1; % time step
final_t = 50;
for t=dt:dt:final_t
    [M_trans, RHS_trans] = transientTerm(c_old, dt, alfa);
    Mconv = convectionTvdTerm(u, c, FL);
    M = M_trans-Mdiff+Mbc+Mconv;
    RHS = RHS_trans+RHSbc;
    c = solvePDE(m,M, RHS);
    c_old = c;
    figure(1);visualizeCells(c);drawnow;
end
%% visualization
 figure(1);visualizeCells(c);
