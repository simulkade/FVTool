%% Transient diffusion equation
% a tutorial adapted from the fipy 1D diffusion example
% see: <http://www.ctcms.nist.gov/fipy/examples/diffusion/index.html>
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
m = createMesh1D(Nx, L);
%% Create the boundary condition structure
BC = createBC(m); % all Neumann boundary condition structure
BC.left.a = 0; BC.left.b=1; BC.left.c=1; % left boundary
BC.right.a = 0; BC.right.b=1; BC.right.c=0; % right boundary
x = m.cellcenters.x;
%% define the transfer coeffs
D_val = 1;
D = createCellVariable(m, D_val);
alfa = createCellVariable(m, 1);
%% define initial values
c_init = 0;
c_old = createCellVariable(m, c_init, BC); % initial values
c = c_old; % assign the old value of the cells to the current values
%% loop
dt = 0.1; % time step
final_t = 100;
for t=dt:dt:final_t
    [M_trans, RHS_trans] = transientTerm(c, dt, alfa);
    Dave = harmonicMean(D);
    Mdiff = diffusionTerm(Dave);
    [Mbc, RHSbc] = boundaryCondition(BC);
    M = M_trans-Mdiff+Mbc;
    RHS = RHS_trans+RHSbc;
    c = solvePDE(m,M, RHS);
    c_analytical = 1-erf(x/(2*sqrt(D_val*t)));
    c_old = c;
%     disp(t)
end
%% visualization
figure(1);plot(x, c.value(2:Nx+1), 'o', x, c_analytical);
xlabel('Length [m]'); ylabel('c');
legend('Numerical', 'Analytical');
