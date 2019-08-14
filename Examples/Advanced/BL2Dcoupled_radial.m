% Coupled nonlinear PDE's
% Buckley Leverett equation
% dependent variables: pressure and water saturation
% Prepared for educational purposes by ** AAE **
clc; clear;
%% define the geometry
well_radius = 0.1; % [m]
x_well_radius = 10*well_radius; % area close to the well
Nx_mid = 50; % number of cells in x direction
Ny_mid = 50; % number of cells in y direction
Lx = 50; % [m]
Ly = 50; % [m]
x = [0:well_radius:x_well_radius linspace(x_well_radius+well_radius, Lx-x_well_radius, Nx_mid) Lx-x_well_radius+well_radius:well_radius:Lx]; 
y = [0:well_radius:x_well_radius linspace(x_well_radius+well_radius, Ly-x_well_radius, Ny_mid) Ly-x_well_radius+well_radius:well_radius:Ly];
m = createMesh2D(x, y); % creates a 2D mesh
Nx = length(x)-1;
Ny = length(y)-1;
%% define the physical parametrs
krw0 = 1.0;
kro0 = 0.76;
nw = 2.4;
no = 2.0;
sor=0.12;
swc=0.09;
sws=@(sw)((sw>swc).*(sw<1-sor).*(sw-swc)/(1-sor-swc)+(sw>=1-sor).*ones(size(sw)));
kro=@(sw)((sw>=swc).*kro0.*(1-sws(sw)).^no+(sw<swc).*(1+(kro0-1)/swc*sw));
krw=@(sw)((sw<=1-sor).*krw0.*sws(sw).^nw+(sw>1-sor).*(-(1-krw0)/sor.*(1.0-sw)+1.0));
dkrwdsw=@(sw)((sw<=1-sor).*nw.*krw0.*(1/(1-sor-swc)).*sws(sw).^(nw-1)+(sw>1-sor)*((1-krw0)/sor));
dkrodsw=@(sw)((sw>=swc).*(-kro0*no*(1-sws(sw)).^(no-1))/(-swc-sor+1)+(sw<swc).*((kro0-1)/swc));
p0 = 100e5; % [bar] pressure
pin = 150e5; % [bar] injection pressure at the left boundary
u_in= 1.0/(24*3600); % [m/s] equal to 1 m/day
sw0 = swc+0.1; % initial water saturation
sw_in = 1;
mu_oil = 2e-3; % [Pa.s] oil viscosity
mu_water = 1e-3; % [Pa.s] water viscosity
% reservoir
k0 = 0.1e-12; % [m^2] average reservoir permeability
phi0 = 0.2; % average porosity
clx=1.2;
cly=0.2;
V_dp=0.1; % Dykstra-Parsons coef.
perm_val= field2d(Nx,Ny,k0,V_dp,clx,cly);
k=createCellVariable(m, perm_val);
phi=createCellVariable(m, phi0);
lw = geometricMean(k)/mu_water;
lo = geometricMean(k)/mu_oil;
Lw = @(sw)(krw(sw));
Lo = @(sw)(k/mu_oil*kro(sw));
dLwdsw = @(sw)(k/mu_water*dkrwdsw(sw));
dLodsw = @(sw)(k/mu_oil*dkrodsw(sw));
%% Define the boundaries
BCp = createBC(m); % Neumann BC for pressure
BCs = createBC(m); % Neumann BC for saturation
% bottom left boundary pressure gradient
BCp.left.a(1)=(krw(sw_in)*lw.xvalue(1,1)+kro(sw_in)*lo.xvalue(1,1)); BCp.left.b(1)=0; BCp.left.c(1)=-u_in;
BCp.bottom.a(1)=(krw(sw_in)*lw.yvalue(1,1)+kro(sw_in)*lo.yvalue(1,1)); BCp.bottom.b(1)=0; BCp.bottom.c(1)=-u_in;
% change the top right boandary to constant pressure (Dirichlet)
BCp.right.a(end)=0; BCp.right.b(end)=1; BCp.right.c(end)=p0;
BCp.top.a(end)=0; BCp.top.b(end)=1; BCp.top.c(end)=p0;
% change the bottom left boundary to constant saturation (Dirichlet)
BCs.left.a(1)=0; BCs.left.b(1)=1; BCs.left.c(1)=1;
BCs.bottom.a(1)=0; BCs.bottom.b(1)=1; BCs.bottom.c(1)=1;
%% define the time step and solver properties
% dt = 1000; % [s] time step
dt=(Lx/Nx)/u_in/10; % [s]

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
dp_alwd= 100.0; % Pa
dsw_alwd= 0.1;
t = 0;
oil_init = domainInt(1-sw_old);
R(1) = 0;
i=0;
v_domain = domainInt(phi);
n_pv = 4; % injecting 4 pv of water
t_end = n_pv*v_domain/(u_in*well_radius*2); % final simulation time
while (t<t_end)
    error_p = 1e5;
    error_sw = 1e5;
    % Implicit loop
    loop_count=0;
    while ((error_p>eps_p) || (error_sw>eps_sw))
        loop_count=loop_count+1;
        if loop_count>10
            break
        end
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
        p_new = reshapeCell(m,full(x(1:(Nx+2)*(Ny+2))));
        sw_new = reshapeCell(m,full(x((Nx+2)*(Ny+2)+1:end)));
        % calculate error values
        error_p = max(abs((p_new(:)-p.value(:))./p_new(:)));
        error_sw = max(abs(sw_new(:)-sw.value(:)));
        % assign new values of p and sw
        p.value = p_new;
        sw.value = sw_new;
    end
    if loop_count>10
      p=p_old;
      sw=sw_old;
      dt=dt/5;
      continue
    end
    dsw=max(abs(sw_new(:)-sw_old.value(:))./sw_new(:));
    t=t+dt
    dt=min([dt*(dsw_alwd/dsw), 2*dt, t_end-t]);
    p_old = p;
    sw_old = sw;
    oil = domainInt(1-sw);
    i=i+1;
    R(i) = (oil_init-oil)/oil_init;
    p_inj(i) = p.value(2,2);
    t_series(i) = t;
%     figure(1);visualizeCells(sw); drawnow;
end
