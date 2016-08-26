% Coupled nonlinear PDE's
% Buckley Leverett equation
% dependent variables: pressure and water saturation
% Prepared for educational purposes by ** AAE **
clc; clear;
%% define the geometry
Nx = 50; % number of cells in x direction
Ny = 30; % number of cells in y direction
W = 50; % [m] length of the domain in x direction
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
q_in=5e-5; % [m^3/s] water injection
sw0 = 0; % initial water saturation
swin = 1;
mu_oil = 1e-3; % [Pa.s] oil viscosity
mu_water = 1e-3; % [Pa.s] water viscosity
% reservoir
k0 = 2e-12; % [m^2] average reservoir permeability
k_frac=100e-12; % [m^2] fracture permeability
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
krw0 = 0.4;
kro0 = 0.8;
nw = 1;
no = 1;
krw = @(sw)(krw0*sw.^nw);
dkrwdsw = @(sw)(krw0*nw*sw.^(nw-1));
kro = @(sw)(kro0*(1-sw).^no);
dkrodsw = @(sw)(-kro0*no*(1-sw).^(no-1));
lw = geometricMean(k)/mu_water;
lo = geometricMean(k)/mu_oil;
Lw = @(sw)(krw(sw));
Lo = @(sw)(k/mu_oil*kro(sw));
dLwdsw = @(sw)(k/mu_water*dkrwdsw(sw));
dLodsw = @(sw)(k/mu_oil*dkrodsw(sw));
%% Define the boundaries
BCp = createBC(m); % Neumann BC for pressure
BCs = createBC(m); % Neumann BC for saturation
% change the right boandary to constant pressure (Dirichlet)
BCp.right.a(end-5:end)=0; BCp.right.b(end-5:end)=1; BCp.right.c(end-5:end)=p0;
% change the left boundary to constant flow
BCp.left.a(1:5)=-k_frac/mu_water; BCp.left.b(1:5)=0; BCp.left.c(1:5)=q_in;
% change the left boundary to constant saturation (Dirichlet)
BCs.left.a(1:5)=0; BCs.left.b(1:5)=1; BCs.left.c(1:5)=1;
%% define the time step and solver properties
dt = 10; % [s] time step
t_end = 1000*dt; % [s] final time
eps_p = 1e-5; % pressure accuracy
eps_sw = 1e-5; % saturation accuracy
%% define the variables
sw_old = createCellVariable(m, sw0, BCs);
p_old = createCellVariable(m, p0, BCp);
sw = sw_old;
p = p_old;
uw = -gradientTerm(p_old); % an estimation of the water velocity
%% start the main loop
% generate intial pressure profile (necessary to initialize the fully
% implicit solver)

t = 0;
while (t<t_end)
    error_p = 1e5;
    error_sw = 1e5;
    % estimation loop: sequential
%     for i = 1:2
%         pgrad = gradientTerm(m, p.value);
%         sw_face = upwindMean(m, -pgrad, sw.value); % average value of water saturation
%         labdao = lo.*funceval(kro, sw_face);
%         labdaw = lw.*funceval(krw, sw_face);
%         labda = labdao+labdaw;
%         Mdiffp1 = diffusionTerm(m, -labda);
%         [Mbcp, RHSbcp] = boundaryCondition(m, BCp);
%         p.value = reshapeCell(m, full((Mdiffp1+Mbcp)\RHSbcp));
%         for j = 1:2
%             pgrad = gradientTerm(m, p.value);
%             sw_face = upwindMean(m, -pgrad, sw.value);
%             labdaw = lw.*funceval(krw, sw_face);
%             dlabdawdsw = lw.*funceval(dkrwdsw, sw_face);
%             [Mbcsw, RHSbcsw] = boundaryCondition(m, BCs);
%             [Mtranssw2, RHStrans2] = transientTerm(m, phi, dt, sw);
%             Mconvsw2 = convectionUpwindTerm(m, -dlabdawdsw.*pgrad);
%             RHSs = RHSbcsw+RHStrans2+divergenceTerm(m, (labdaw-dlabdawdsw.*sw_face).*pgrad);
%             Ms = Mbcsw+Mconvsw2+Mtranssw2;
%             sw.value = reshapeCell(m, full(Ms\RHSs));
%         end
%     end
    % Implicit loop
    hold all
    while ((error_p>eps_p) || (error_sw>eps_sw))
        % calculate parameters
        pgrad = gradientTerm(p);
        sw_face = upwindMean(sw, -pgrad); % average value of water saturation
        labdao = lo.*funceval(kro, sw_face);
        labdaw = lw.*funceval(krw, sw_face);
        dlabdaodsw = lo.*funceval(dkrodsw, sw_face);
        dlabdawdsw = lw.*funceval(dkrwdsw, sw_face);
        labda = labdao+labdaw;
        dlabdadsw = dlabdaodsw+dlabdawdsw;
        % compute [Jacobian] matrices
        Mdiffp1 = diffusionTerm(-labda);
        Mdiffp2 = diffusionTerm(-labdaw);
        Mconvsw1 = convectionUpwindTerm(-dlabdadsw.*pgrad);
        Mconvsw2 = convectionUpwindTerm(-dlabdawdsw.*pgrad);
        [Mtranssw2, RHStrans2] = transientTerm(sw_old, dt, phi);
        % Compute RHS values
        RHS1 = divergenceTerm(-dlabdadsw.*sw_face.*pgrad);
        RHS2 = divergenceTerm(-dlabdawdsw.*sw_face.*pgrad);
        % include boundary conditions
        [Mbcp, RHSbcp] = boundaryCondition(BCp);
        [Mbcsw, RHSbcsw] = boundaryCondition(BCs);
        % Couple the equations; BC goes into the block on the main diagonal
        M = [Mdiffp1+Mbcp Mconvsw1; Mdiffp2 Mconvsw2+Mtranssw2+Mbcsw];
        RHS = [RHS1+RHSbcp; RHS2+RHStrans2+RHSbcsw];
        % solve the linear system of equations
        x = M\RHS;
        % x = agmg(M, RHS, [], 1e-10, 500, [], [p.value(:); sw.value(:)]);
        % separate the variables from the solution
        p_new = reshapeCell(m,full(x(1:(2*Nx+1)*(2*Ny+1))));
        sw_new = reshapeCell(m,full(x((2*Nx+1)*(2*Ny+1)+1:end)));
        % calculate error values
        error_p = max(abs((p_new(:)-p.value(:))./p_new(:)));
        error_sw = max(abs(sw_new(:)-sw.value(:)));
        % assign new values of p and sw
        p.value = p_new;
        sw.value = sw_new;
    end
    t=t+dt;
    p_old = p;
    sw_old = sw;
    figure(1);visualizeCells(sw);shading interp; drawnow;
    sw_tot=sw.value.*X.*Y/(W*H);
    oil_prod=(sum(sum(sw_tot(2:end-1,2:end-1))));
    figure(2);semilogx(t, p_new(2,2), 'o'); drawnow;
end
