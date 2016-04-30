% Coupled nonlinear PDE's
% Buckley Leverett equation
% dependent variables: pressure and water saturation
% Prepared for educational purposes by ** AAE **
clc; clear;
%% define the geometry
Nx = 50; % number of cells in x direction
Ny = 50; % number of cells in y direction
Nxin=Nx-19;
Nyin=Ny-19;
W = 0.5; % [m] length of the domain in x direction
H = 0.5; % [m] length of the domain in y direction
% m = createMesh1D(Nx, W);
m = createMesh2D(Nx, Ny, W, H); % creates a 2D mesh
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
sw0=ones(Nx, Ny);
sw0(10:end-10, 10:end-10)=swc+0.2;
% sw0 = swc+0.1; % initial water saturation
sw_in = 1;
mu_oil = 2e-3; % [Pa.s] oil viscosity
mu_water = 1e-3; % [Pa.s] water viscosity
% reservoir
k0 = 0.1e-12; % [m^2] average reservoir permeability
k0_frac=100e-12; % [m^2]
phi0 = 0.2; % average porosity
phi0_frac=1.0;
teta_ow=deg2rad(30);
gama_ow=0.03; % N/m
labda=10.0;
eps1=1e-6;
clx=0.05;
cly=0.05;
V_dp=0.1; % Dykstra-Parsons coef.
perm_valin= k0;%field2d(Nxin,Nyin,k0,V_dp,clx,cly);
perm_val=k0_frac*ones(Nx, Ny);
perm_val(10:end-10, 10:end-10)=perm_valin;
k=createCellVariable(m, perm_val);
phi_val=phi0_frac*ones(Nx, Ny);
phi_val(10:end-10, 10:end-10)=phi0;
phi=createCellVariable(m, phi_val);
pce=gama_ow*cos(teta_ow)*(phi./k).^0.5;
pce_face=gama_ow*cos(teta_ow)*(arithmeticMean(phi)./geometricMean(k)).^0.5;
sco=swc+0.01;
pc=@(sw)(pce.*(sws(sw)+eps).^(-1/labda)); % it can also be defined for each block
dpc=@(sw)((-1/labda)*(1/(1-sor-swc)).*pce_face.*(sws(sw)+eps).^(-1/labda-1));
dpcdk=@(sw)(0.5*gama_ow*cos(teta_ow)*(geometricMean(k)./arithmeticMean(phi)).^0.5.*(sws(sw)+eps).^(-1/labda));
grad_phik=gradientTerm(phi./k);
% sw_plot=linspace(0,1, 10000);
% plot(sw_plot, pc(sw_plot))
lw = geometricMean(k)/mu_water;
lo = geometricMean(k)/mu_oil;
Lw = @(sw)(krw(sw));
Lo = @(sw)(k/mu_oil*kro(sw));
dLwdsw = @(sw)(k/mu_water*dkrwdsw(sw));
dLodsw = @(sw)(k/mu_oil*dkrodsw(sw));
%% Define the boundaries: all fixed Sw=1, fixed pressure everywhere(?)
BCp = createBC(m); % Neumann BC for pressure
BCs = createBC(m); % Neumann BC for saturation
% left boundary pressure gradient
% BCp.left.a(:)=(krw(sw_in)*lw.xvalue(1,:)+kro(sw_in)*lo.xvalue(1,:)); BCp.left.b(:)=0; BCp.left.c(:)=-u_in;
% change the right boandary to constant pressure (Dirichlet)
% BCp.left.a(:)=0; BCp.left.b(:)=1; BCp.left.c(:)=p0;
% BCp.right.a(:)=0; BCp.right.b(:)=1; BCp.right.c(:)=p0;
BCp.top.a(:)=0; BCp.top.b(:)=1; BCp.top.c(:)=p0;
% BCp.bottom.a(:)=0; BCp.bottom.b(:)=1; BCp.bottom.c(:)=p0+100;
% change the left boundary to constant saturation (Dirichlet)
% BCs.left.a(:)=0; BCs.left.b(:)=1; BCs.left.c(:)=1.0-sor;
% BCs.right.a(:)=0; BCs.right.b(:)=1; BCs.right.c(:)=1.0-sor;
% BCs.top.a(:)=0; BCs.top.b(:)=1; BCs.top.c(:)=1.0;
BCs.bottom.a(:)=0; BCs.bottom.b(:)=1; BCs.bottom.c(:)=1.0;
%% define the time step and solver properties
% dt = 1000; % [s] time step
% dt=(W/Nx)/u_in/20; % [s]
dt=0.3;
t_end = 1000*dt; % [s] final time
eps_p = 1e-7; % pressure accuracy
eps_sw = 1e-7; % saturation accuracy
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
% for i=1:5
    error_p = 1e5;
    error_sw = 1e5;
    % Implicit loop
%     while ((error_p>eps_p) || (error_sw>eps_sw))
    for j=1:3
        % calculate parameters
        pgrad = gradientTerm(p);
%         pcgrad=gradientTerm(pc(sw));
        sw_face = upwindMean(sw, -pgrad); % average value of water saturation
        sw_grad=gradientTerm(sw);
        sw_ave=arithmeticMean(sw);
        pcgrad=dpc(sw_ave).*sw_grad+dpcdk(sw_ave).*grad_phik;
        % solve for pressure at known Sw
        labdao = lo.*funceval(kro, sw_face);
        labdaw = lw.*funceval(krw, sw_face);
%         dlabdaodsw = lo.*funceval(dkrodsw, sw_face);
%         dlabdawdsw = lw.*funceval(dkrwdsw, sw_face);
        labda = labdao+labdaw;
        % compute [Jacobian] matrices
        Mdiffp1 = diffusionTerm(-labda);
        RHSpc1=divergenceTerm(labdao.*pcgrad);
        [Mbcp, RHSbcp] = boundaryCondition(BCp);
        RHS1 = RHSpc1+RHSbcp; % with capillary
        p_new=solvePDE(m, Mdiffp1+Mbcp, RHS1);
        
        % solve for Sw
        pgrad = gradientTerm(p_new);
        uw=-labdaw.*pgrad;
        [Mbcsw, RHSbcsw] = boundaryCondition(BCs);
        RHS_sw=-divergenceTerm(uw);
        sw_new=solveExplicitPDE(sw_old, dt, RHS_sw, BCs, phi);

        error_p = max(abs((p_new.value(:)-p.value(:))./p_new.value(:)))
        error_sw = max(abs(sw_new.value(:)-sw.value(:)))
        % assign new values of p and sw
        p = p_new;
        sw = sw_new;
    end
    t=t+dt;
    p_old = p;
    sw_old = sw;
    figure(1);visualizeCells(1-sw); drawnow;
end