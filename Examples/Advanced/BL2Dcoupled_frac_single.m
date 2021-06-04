% single phase flow in fractured reservoir
% dependent variables: pressure
% Prepared for educational purposes by Ali Akbar Eftekhari
% DTU, a cold Friday afternoon
% Does not work well
% Written by Ali A. Eftekhari
% Last checked: June 2021
clc; clear;
%% define the geometry
Nx = 50; % number of cells in x direction
Ny = 10; % number of cells in y direction
W = 500; % [m] length of the domain in x direction
H = 30; % [m] length of the domain in y direction
x1=linspace(0,W, Nx);
x2=x1+0.001;
x=zeros(2*Nx,1);
j=0;
for i=1:Nx
    j=j+1;
    x(j)=x1(i);
    j=j+1;
    x(j)=x2(i);
end
y1=linspace(0,H, Ny);
y2=y1+0.01;
y=zeros(2*Ny,1);
j=0;
for i=1:Ny
    j=j+1;
    y(j)=y1(i);
    j=j+1;
    y(j)=y2(i);
end
m = createMesh2D(x, y); % creates a 2D mesh
[X,Y]=ndgrid(m.cellsize.x, m.cellsize.y);
%% define the physical parametrs
p0 = 1e5; % [bar] pressure
pin = 50e5; % [bar] injection pressure at the left boundary
q_in=1e-6; % [m^3/s] water injection
mu_water = 1e-3; % [Pa.s] water viscosity
c1=1e-3; % compressibility
% reservoir
k0 = 1e-12; % [m^2] average reservoir permeability
k_frac=1000e-12; % [m^2] fracture permeability
phi0 = 0.2; % average porosity
clx=0.05;
cly=0.05;
V_dp=0.1; % Dykstra-Parsons coef.
perm_val= field2d(2*Nx,2*Ny,k0,V_dp,clx,cly);
k=createCellVariable(m, k0);
for i=4:2:2*Nx
    k.value(i,3:end-2)=k_frac;
end
for j=4:2:2*Ny
    k.value(3:end-2,j)=k_frac;
end
k.value(1:2,:)=k_frac;
k.value(end-1:end,:)=k_frac;
k.value(:,1:2)=k_frac;
k.value(:,end-1:end)=k_frac;
%k.value(1:2,:)=k_frac;
phi=createCellVariable(m, phi0);
k_face=harmonicMean(k);
labda_face=k_face/mu_water;
Mdiffp = diffusionTerm(-labda_face);

%% Define the boundaries
BCp = createBC(m); % Neumann BC for pressure
% change the right boandary to constant pressure (Dirichlet)
BCp.right.a(end-5:end)=0; BCp.right.b(end-5:end)=1; BCp.right.c(end-5:end)=p0;
% change the left boundary to constant flow
BCp.left.a(1:5)=-k_frac/mu_water; BCp.left.b(1:5)=0; BCp.left.c(1:5)=q_in;
[Mbcp, RHSbcp] = boundaryCondition(BCp);
%% define the time step and solver properties
dt = 1000; % [s] time step
t_end = 1000*dt; % [s] final time
%% define the variables
p_old = createCellVariable(m, p0, BCp);
p = p_old;
uw = -gradientTerm(p_old); % an estimation of the water velocity
%% start the main loop
% generate intial pressure profile (necessary to initialize the fully
% implicit solver)

t = 0;
while (t<t_end)
    hold all
    % compute [Jacobian] matrices
    [Mtrans, RHStrans] = transientTerm(p_old, dt, c1*phi);
    M = Mdiffp+Mbcp+Mtrans;
    RHS = RHSbcp+RHStrans;
    % solve the linear system of equations
    p_new = solvePDE(m, M, RHS);
    t=t+dt;
    p_old = p_new;
    figure(1);visualizeCells(p_new);shading interp; drawnow;
    figure(2);semilogx(t, p_new.value(2,2), 'o'); drawnow;
end
