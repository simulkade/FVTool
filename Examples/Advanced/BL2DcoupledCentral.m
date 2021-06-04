% Coupled nonlinear PDE's
% Buckley Leverett equation for variables pressure and water saturation
% Prepared for educational purposes by ** AAE **
% using a second order scheme to make sharper fronts; quickly fails
% Written by Ali A. Eftekhari
% Last checked: June 2021
clc
%% define the geometry
Nx = 150; % number of cells in x direction
Ny = 100; % number of cells in y direction
W = 50; % [m] length of the domain in x direction
H = 30; % [m] length of the domain in y direction
m = createMesh2D(Nx, Ny, W, H); % creates a 2D mesh
%% define the physical parametrs
p0 = 1e5; % [bar] pressure
pin = 50e5; % [bar] injection pressure at the left boundary
sw0 = 0; % initial water saturation
swin = 1;
mu_oil = 1; % [Pa.s] oil viscosity
mu_water = 1e-3; % [Pa.s] water viscosity
% reservoir
k0 = 1e-12; % [m^2] average reservoir permeability
phi0 = 0.2; % average porosity
clx=0.05;
cly=0.05;
V_dp=0.6; % Dykstra-Parsons coef.
perm_val= field2d(Nx,Ny,k0,V_dp,clx,cly);
k=createCellVariable(m, perm_val);
phi=createCellVariable(m, phi0);
krw0 = 1;
kro0 = 1;
nw = 4;
no = 2;
krw = @(sw)(krw0*sw.^nw);
dkrwdsw = @(sw)(krw0*nw*sw.^(nw-1));
kro = @(sw)(kro0*(1-sw).^no);
dkrodsw = @(sw)(-kro0*no*(1-sw).^(no-1));
lw = harmonicMean(k)/mu_water;
lo = harmonicMean(k)/mu_oil;
Lw = @(sw)(krw(sw));
Lo = @(sw)(k/mu_oil*kro(sw));
dLwdsw = @(sw)(k/mu_water*dkrwdsw(sw));
dLodsw = @(sw)(k/mu_oil*dkrodsw(sw));
%% Define the boundaries
BCp = createBC(m); % Neumann BC for pressure
BCs = createBC(m); % Neumann BC for saturation
% change the left and right boandary to constant pressure (Dirichlet)
BCp.left.a(:)=0; BCp.left.b(:)=1; BCp.left.c(:)=pin;
BCp.right.a(:)=0; BCp.right.b(:)=1; BCp.right.c(:)=p0;
% change the left boundary to constant saturation (Dirichlet)
BCs.left.a(:)=0; BCs.left.b(:)=1; BCs.left.c(:)=1;
%% define the time step and solver properties
dt = 50000; % [s] time step
t_end = 10000*dt; % [s] final time
eps_p = 1e-7; % pressure accuracy
eps_sw = 1e-10; % saturation accuracy
%% define the variables
sw_old = createCellVariable(m, sw0, BCs);
p_old = createCellVariable(m, p0, BCp);
sw = sw_old;
p = p_old;
uw = -gradientTerm(p); % an estimation of the water velocity
%% start the main loop
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
    while ((error_p>eps_p) || (error_sw>eps_sw))
        % calculate parameters
        pgrad = gradientTerm(p);
        sw_face = arithmeticMean(sw); % average value of water saturation
        labdao = lo.*funceval(kro, sw_face);
        labdaw = lw.*funceval(krw, sw_face);
        dlabdaodsw = lo.*funceval(dkrodsw, sw_face);
        dlabdawdsw = lw.*funceval(dkrwdsw, sw_face);
        labda = labdao+labdaw;
        dlabdadsw = dlabdaodsw+dlabdawdsw;
        % compute [Jacobian] matrices
        Mdiffp1 = diffusionTerm(-labda);
        Mdiffp2 = diffusionTerm(-labdaw);
        Mconvsw1 = convectionTerm(-dlabdadsw.*pgrad);
        Mconvsw2 = convectionTerm(-dlabdawdsw.*pgrad);
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
    t=t+dt;
    p_old = p;
    sw_old = sw;
    figure(1);visualizeCells(sw);drawnow;
end
