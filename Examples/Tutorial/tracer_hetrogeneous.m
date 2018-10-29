% see http://fvt.simulkade.com/posts/2015-05-11-tracer-flow-porous-media.html
Nx=50;
Ny=50;
Lx=1.0; % (m)
Ly=1.0; % (m)
% physical values
D_val=1.0e-9; % (m^2/s)
mu_val=1e-3; % (Pa.s)
poros=0.2;
perm_val=1.0e-12; % (m^2)
clx=0.05;
cly=0.05;
V_dp=0.6;
perm= field2d(Nx,Ny,perm_val,V_dp,clx,cly);
% physical system
u_inj=1.0/(3600*24); % (m/s)
c_init=0.0;
c_inj=1.0;
p_out=100e5; % (Pa)
% create mesh and assign values to the domain
m= createMesh2D(Nx, Ny, Lx, Ly);
k=createCellVariable(m, perm);
phi=createCellVariable(m,poros);
D=createCellVariable(m, D_val);
% Define the boundaries
BCp = createBC(m); % Neumann BC for pressure
BCc = createBC(m); % Neumann BC for concentration
% change the right boandary to constant pressure (Dirichlet)
BCp.right.a(:)=0.0;
BCp.right.b(:)=1.0;
BCp.right.c(:)=p_out;
% left boundary
BCp.left.a(:)=-perm_val/mu_val;
BCp.left.c(:)=u_inj;
% change the left boundary to constant concentration (Dirichlet)
BCc.left.a(:)=0.0;
BCc.left.b(:)=1.0;
BCc.left.c(:)=1.0;
labda_face=harmonicMean(k/mu_val);
Mdiffp=diffusionTerm(labda_face);
[Mbcp, RHSp] = boundaryCondition(BCp);
Mp= Mdiffp+Mbcp;
p_val=solvePDE(m, Mp, RHSp);
figure(1)
visualizeCells(p_val);
title('Pressure profile (Pa)');
colorbar();

% estimate a practical time step
n_loop=50;
dt=Lx/u_inj/(5*100); % (s)
% find the velocity vector
u=-labda_face.*gradientTerm(p_val); % (m/s)
% find the matrices of coefficients
D_face=harmonicMean(phi.*D);
Mdiffc=diffusionTerm(D_face);
Mconvc=convectionUpwindTerm(u);
[Mbcc, RHSbcc] = boundaryCondition(BCc);
% initialize
c_old = createCellVariable(m, c_init, BCc);
c_val = c_old;
% start the loop
for i=1:n_loop
    [Mtransc, RHStransc] = transientTerm(c_old, dt, phi);
    Mc=-Mdiffc+Mconvc+Mbcc+Mtransc;
    RHSc=RHSbcc+RHStransc;
    c_val=solvePDE(m, Mc, RHSc);
    c_old=c_val;
end
figure(2);
subplot(1,2,1);
visualizeCells(c_val);
title('concentration profile');
colorbar();
subplot(1,2,2);
visualizeCells(k);
colorbar()
title('perm field');
