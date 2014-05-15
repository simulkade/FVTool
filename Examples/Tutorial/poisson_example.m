%% Poisson equation
% see this link
% http://scicomp.stackexchange.com/questions/8577/peculiar-error-when-solving-the-poisson-equation-on-a-non-uniform-mesh-1d-only
% Strange behavior when change the number of grids from even to odd
% Wrong results does not always mean that the code has bugs. 
% Wrong use of the code can also give you wrong results
clc; clear;

%% Define the domain and create a mesh structure
L = 20;  % domain length
Nx = 10000; % number of cells
m = buildMesh1D(Nx, L);
%% Create the boundary condition structure
BC = createBC(m); % all Neumann boundary condition structure
BC.left.a = 0; BC.left.b=1; BC.left.c=0; % left boundary
% BC.right.a = 0; BC.right.b=1; BC.right.c=0; % right boundary
x = m.cellcenters.x-10; % move the domain to [-10,10]
%% define the transfer coeffs
D_val = 1;
D = createFaceVariable(m, D_val);
%% define source term
rho = @(x)(-1.0*((x>=-1.0)&(x<=0))+((x>0)&(x<=1)));
s1 = sourceExplicitTerm(m, rho(x));
Mdiff = diffusionTerm(m, D);
[Mbc, RHSbc] = boundaryCondition(m, BC);
M = Mdiff+Mbc;
RHS = -s1+RHSbc;
c = solvePDE(m,M, RHS);
%% visualization
figure(1);plot(x, c(2:Nx+1), x, rho(x));
xlabel('Length [m]'); ylabel('c');
legend('Numerical', 'charge');
