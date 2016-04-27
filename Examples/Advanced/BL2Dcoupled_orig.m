% Coupled nonlinear PDE's
% Buckley Leverett equation
% dependent variables: pressure and water saturation
% Prepared for educational purposes by ** AAE **
clc; clear;
%% define the geometry
Nx = 100; % number of cells in x direction
Ny = 30; % number of cells in y direction
W = 300; % [m] length of the domain in x direction
H = 30; % [m] length of the domain in y direction
m = createMesh1D(Nx, W);
% m = createMesh2D(Nx, Ny, W, H); % creates a 2D mesh
%% define the physical parametrs
krw0_ww = 0.3;
krw0_ow = 1.0;
kro0_ww = 0.6;
kro0_ow = 0.76;
nw = 2.4;
no = 2.0;
sor_ww=0.12;
sor_ow=0.12;
swc_ww=0.09;
swc_ow=0.09;
SF=createFaceVariable(m, 1.0); % 1 is water wet, 0 is oil wet
krw0=krw0_ww*SF+krw0_ow*(1-SF);
kro0=kro0_ww*SF+kro0_ow*(1-SF);
sor=sor_ww*SF+sor_ow*(1-SF);
swc=swc_ww*SF+swc_ow*(1-SF);
sws=@(sw, sor, swc)((sw>swc).*(sw<1-sor).*(sw-swc)./(1-sor-swc)+(sw>=1-sor));
kro=@(sw, kro0, sor, swc)((sw>=swc).*kro0.*(1-sws(sw, sor)).^no+(sw<swc).*(1+(kro0-1)/swc*sw));
krw=@(sw, krw0, sor, swc)((sw<=1-sor).*krw0.*sws(sw, sor).^nw+(sw>1-sor).*(-(1-krw0)/sor.*(1.0-sw)+1.0));
dkrwdsw=@(sw, krw0, sor, swc)((sw<=1-sor).*nw.*krw0.*(1/(1-sor-swc)).*sws(sw, sor).^(nw-1)+(sw>1-sor).*((1-krw0)./sor));
dkrodsw=@(sw, kro0, sor, swc)((sw>=swc).*(-kro0*no*(1-sws(sw, sor)).^(no-1))/(-swc-sor+1)+(sw<swc).*((kro0-1)./swc));
p0 = 100e5; % [bar] pressure
pin = 150e5; % [bar] injection pressure at the left boundary
u_in= 1.0/(24*3600); % [m/s] equal to 1 m/day
sw0 = swc+0.1; % initial water saturation
sw_in = 1;
mu_oil = 2e-3; % [Pa.s] oil viscosity
mu_water = 1e-3; % [Pa.s] water viscosity
% reservoir
k0 = 2e-12; % [m^2] average reservoir permeability
phi0 = 0.2; % average porosity
teta_ow=deg2rad(30);
gama_ow=0.03; % N/m
labda=10.0;
eps1=1e-7;
clx=0.2;
cly=0.2;
V_dp=0.7; % Dykstra-Parsons coef.
if m.dimension<2 % 1D model
    perm_val=k0;
elseif m.dimension<3 % 2D model
    perm_val= field2d(Nx,Ny,k0,V_dp,clx,cly);
end
k=createCellVariable(m, perm_val);
phi=createCellVariable(m, phi0);
n_cells=numel(phi.value);
pce=gama_ow*cos(teta_ow)*(phi./k).^0.5;
pce_face=gama_ow*cos(teta_ow)*(arithmeticMean(phi)./geometricMean(k)).^0.5;
pc=@(sw)(pce.*(sws(sw)+eps1).^(-1/labda)); % it can also be defined for each block
dpc=@(sw)((-1/labda)*(1/(1-sor-swc)).*pce_face.*(sws(sw)+eps1).^(-1/labda-1));
dpcdk=@(sw)(0.5*gama_ow*cos(teta_ow)*(geometricMean(k)./arithmeticMean(phi)).^0.5.*(sws(sw)+eps1).^(-1/labda));
grad_phik=gradientTerm(phi./k);
% sw_plot=linspace(0,1, 10000);
% plot(sw_plot, pc(sw_plot))
lw = geometricMean(k)/mu_water;
lo = geometricMean(k)/mu_oil;
Lw = @(sw)(krw(sw));
Lo = @(sw)(k/mu_oil*kro(sw));
dLwdsw = @(sw)(k/mu_water*dkrwdsw(sw));
dLodsw = @(sw)(k/mu_oil*dkrodsw(sw));
%% Define the boundaries
BCp = createBC(m); % Neumann BC for pressure
BCs = createBC(m); % Neumann BC for saturation
% left boundary pressure gradient
BCp.left.a(:)=(krw(sw_in)*lw.xvalue(1,:)+kro(sw_in)*lo.xvalue(1,:)); BCp.left.b(:)=0; BCp.left.c(:)=-u_in;
% change the right boandary to constant pressure (Dirichlet)
% BCp.left.a(:)=0; BCp.left.b(:)=1; BCp.left.c(:)=pin;
BCp.right.a(:)=0; BCp.right.b(:)=1; BCp.right.c(:)=p0;
% change the left boundary to constant saturation (Dirichlet)
BCs.left.a(:)=0; BCs.left.b(:)=1; BCs.left.c(:)=1;
% specific boundary conditions for capillary pressure
BCpc=createBC(m);
% top and bottom are closed. pc gradient is zero (default)
% on the left side it is also zero (assumption, sw=1)
% on the right side, end effect. pc=0 at the right ghost cell.
% BCpc.right.a(:)=0; BCpc.right.b(:)=1; BCpc.right.c(:)=0;
%% define the time step and solver properties
% dt = 1000; % [s] time step
dt=(W/Nx)/u_in/20; % [s]
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
    % Implicit loop
    while ((error_p>eps_p) || (error_sw>eps_sw))
        % calculate parameters
        pgrad = gradientTerm(p);
        sw_face = upwindMean(sw, -pgrad); % average value of water saturation
        sw_grad=gradientTerm(sw);
        sw_ave=arithmeticMean(sw);
        pcgrad=dpc(sw_ave).*sw_grad+dpcdk(sw_ave).*grad_phik;
%         pc_cell= pc(sw);
%         pc_cell.value(end,:)=0; % right ghost cells, end effect.
        % pc(1,:)=0;
%         pcgrad=gradientTerm(pc_cell);
        
        labdao = lo.*funceval(kro, sw_face);
        labdaw = lw.*funceval(krw, sw_face);
        dlabdaodsw = lo.*funceval(dkrodsw, sw_face);
        dlabdawdsw = lw.*funceval(dkrwdsw, sw_face);
        labda = labdao+labdaw;
        dlabdadsw = dlabdaodsw+dlabdawdsw;
        % compute [Jacobian] matrices
        Mdiffp1 = diffusionTerm(-labda);
        Mdiffp2 = diffusionTerm(-labdaw);
%         Mconvsw1 = convectionUpwindTerm(-dlabdadsw.*pgrad); % without capillary
        Mconvsw1 = convectionUpwindTerm(-dlabdawdsw.*pgrad-dlabdaodsw.*(pgrad+pcgrad)); % with capillary
        Mconvsw2 = convectionUpwindTerm(-dlabdawdsw.*pgrad);
        [Mtranssw2, RHStrans2] = transientTerm(sw_old, dt, phi);
        % Compute RHS values
        RHSpc1=divergenceTerm(labdao.*pcgrad);
        RHS1 = divergenceTerm((-dlabdawdsw.*pgrad-dlabdaodsw.*(pgrad+pcgrad)).*sw_face); % with capillary
        RHS2 = divergenceTerm(-dlabdawdsw.*sw_face.*pgrad);
        % include boundary conditions
        [Mbcp, RHSbcp] = boundaryCondition(BCp);
        [Mbcsw, RHSbcsw] = boundaryCondition(BCs);
        % Couple the equations; BC goes into the block on the main diagonal
        M = [Mdiffp1+Mbcp Mconvsw1; Mdiffp2 Mconvsw2+Mtranssw2+Mbcsw];
        RHS = [RHS1+RHSpc1+RHSbcp; RHS2+RHStrans2+RHSbcsw];
        % solve the linear system of equations
        x = M\RHS;
        % x = agmg(M, RHS, [], 1e-10, 500, [], [p.value(:); sw.value(:)]);
        % separate the variables from the solution
        p_new = reshapeCell(m,full(x(1:n_cells)));
        sw_new = reshapeCell(m,full(x(n_cells+1:end)));
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
    figure(1);visualizeCells(1-sw);title([num2str(t/3600/24) ' day']); drawnow;
end
