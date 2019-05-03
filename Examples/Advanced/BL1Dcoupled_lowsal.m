% Coupled nonlinear PDE's
% Buckley Leverett equation
% dependent variables: pressure and water saturation
% Prepared for educational purposes by ** AAE **
clc; clear;
%% define the geometry
Nx = 50; % number of cells in x direction
Ny = 30; % number of cells in y direction
W = 300; % [m] length of the domain in x direction
H = 30; % [m] length of the domain in y direction
m = createMesh1D(Nx, W); % creates a 1D mesh
%% define the initial and boundary conditions for concentrations
c_inj_ion = 1.0; % injection concentration
c_init_ion = 0.0;
BCc = createBC(m);
BCc.left.a = 0.0; BCc.left.b = 1.0; BCc.left.c = c_inj_ion; % left boundary
c_ion_old = createCellVariable(m, c_init_ion, BCc);
c_ion = c_ion_old;
[M_bc_c, RHS_bc_c]=boundaryCondition(BCc);

%% define the diffusion adsorption domains for each cell
% partameters
rho_s = 2700; % kg/m3
a_s = 2000; % m2/kg
a = 1e9; % m^2/m^3
k_lang = 1e-7; % Langmuir adsorption coefficient
betta = 1.0; % Langmuir adsorption coefficient 2
% Define the domain and create a mesh structure
L = 1e-6; % [m]  % domain length
N = 15; % number of cells
m_diff = createMesh1D(N, L);
% Create the boundary condition structure
BC = cell(Nx, 1);
for i = 1:Nx
    BC{i} = createBC(m_diff); % all Neumann boundary condition structure
    BC{i}.left.a = 0; BC{i}.left.b=1; BC{i}.left.c=0; % left boundary
    BC{i}.right.a = 0; BC{i}.right.b=1; BC{i}.right.c=0; % right boundary
end
D_val = 1e-13; % m^2/s effective diffusivity
D = createCellVariable(m_diff, D_val);
D_face = harmonicMean(D);
Mdiff = diffusionTerm(D_face);
alfa = createCellVariable(m_diff, 1);
% initial condition
c_init = 0;
c_old = cell(Nx);
for i = 1:Nx
    c_old{i} = createCellVariable(m_diff, c_init, BC{i}); % initial values
end
c = c_old; % assign the old value of the cells to the current values
%% define the physical parametrs
% define the oil-wet rel-perm
krw0 = 0.5;
kro0 = 0.76;
nw = 2.4;
no = 2.0;
sor=0.15;
swc=0.09;
sws=@(sw)((sw>swc).*(sw<1-sor).*(sw-swc)/(1-sor-swc)+(sw>=1-sor).*ones(size(sw)));
kro=@(sw)((sw>=swc).*kro0.*(1-sws(sw)).^no+(sw<swc).*(1+(kro0-1)/swc*sw));
krw=@(sw)((sw<=1-sor).*krw0.*sws(sw).^nw+(sw>1-sor).*(-(1-krw0)/sor.*(1.0-sw)+1.0));
dkrwdsw=@(sw)((sw<=1-sor).*nw.*krw0.*(1/(1-sor-swc)).*sws(sw).^(nw-1)+(sw>1-sor)*((1-krw0)/sor));
dkrodsw=@(sw)((sw>=swc).*(-kro0*no*(1-sws(sw)).^(no-1))/(-swc-sor+1)+(sw<swc).*((kro0-1)/swc));

% define water wet relperms
krw0_ww = 0.3;
kro0_ww = 0.96;
nw_ww = 2.4;
no_ww = 2.0;
sor_ww =0.05;
swc_ww =0.09;
sws_ww=@(sw)((sw>swc_ww).*(sw<1-sor_ww).*(sw-swc_ww)/(1-sor_ww-swc_ww)+(sw>=1-sor_ww).*ones(size(sw)));
kro_ww=@(sw)((sw>=swc_ww).*kro0_ww.*(1-sws_ww(sw)).^no_ww+(sw<swc_ww).*(1+(kro0_ww-1)/swc_ww*sw));
krw_ww=@(sw)((sw<=1-sor_ww).*krw0_ww.*sws_ww(sw).^nw_ww+(sw>1-sor_ww).*(-(1-krw0_ww)/sor_ww.*(1.0-sw)+1.0));
dkrwdsw_ww=@(sw)((sw<=1-sor_ww).*nw.*krw0_ww.*(1/(1-sor-swc_ww)).*sws_ww(sw).^(nw_ww-1)+(sw>1-sor_ww)*((1-krw0_ww)/sor_ww));
dkrodsw_ww=@(sw)((sw>=swc_ww).*(-kro0_ww*no*(1-sws_ww(sw)).^(no_ww-1))/(-swc_ww-sor_ww+1)+(sw<swc_ww).*((kro0_ww-1)/swc_ww));

%% wettability modifier
SF=createCellVariable(m, 0.0); % 1 is water wet, 0 is oil wet
SF_face = arithmeticMean(SF);


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
clx=1.2;
cly=0.2;
V_dp=0.7; % Dykstra-Parsons coef.
perm_val= k0; %field2d(Nx,Ny,k0,V_dp,clx,cly);
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
% left boundary pressure gradient
BCp.left.a(:)=(krw(sw_in)*lw.xvalue(1,:)+kro(sw_in)*lo.xvalue(1,:)); BCp.left.b(:)=0; BCp.left.c(:)=-u_in;
% change the right boandary to constant pressure (Dirichlet)
% BCp.left.a(:)=0; BCp.left.b(:)=1; BCp.left.c(:)=pin;
BCp.right.a(:)=0; BCp.right.b(:)=1; BCp.right.c(:)=p0;
% change the left boundary to constant saturation (Dirichlet)
BCs.left.a(:)=0; BCs.left.b(:)=1; BCs.left.c(:)=1;
%% define the time step and solver properties
% dt = 1000; % [s] time step
dt=(W/Nx)/u_in/10; % [s]
t_end = 400*dt; % [s] final time
eps_p = 1e-5; % pressure accuracy
eps_sw = 1e-5; % saturation accuracy
%% define the variables
sw_old = createCellVariable(m, sw0, BCs);
oil_init=domainInt(1-sw_old);

p_old = createCellVariable(m, p0, BCp);
sw = sw_old;
p = p_old;
uw = -gradientTerm(p_old); % an estimation of the water velocity
%% start the main loop
% generate intial pressure profile (necessary to initialize the fully
% implicit solver)
rec_fact=0;
t_day=0;
dp_alwd= 100.0; % Pa
dsw_alwd= 0.05;
t = 0;
% fprintf(1, 'progress (%%):  ');
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
        labdao = ((1-SF_face).*funceval(kro, sw_face)+SF_face.*funceval(kro_ww, sw_face)).*lo;
        labdaw = ((1-SF_face).*funceval(krw, sw_face)+SF_face.*funceval(krw_ww, sw_face)).*lw;
        dlabdaodsw = ((1-SF_face).*funceval(dkrodsw, sw_face)+SF_face.*funceval(dkrodsw_ww, sw_face)).*lo;
        dlabdawdsw = ((1-SF_face).*funceval(dkrwdsw, sw_face)+SF_face.*funceval(dkrwdsw_ww, sw_face)).*lw;
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
        p_new = reshapeCell(m,full(x(1:(Nx+2))));
        sw_new = reshapeCell(m,full(x((Nx+2)+1:end)));
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
    % solve the ions flow in the domain
    uw = -labdaw.*pgrad; % water velocity
    Mconv = convectionUpwindTerm(uw);
    dswdt=(sw-sw_old)/dt;
    Msc=linearSourceTerm(phi.*dswdt);
    for j = 1:3
        [M_trans, RHS_trans] = transientTerm(c_ion_old, dt, sw.*phi+a_s*rho_s*k_lang.*(1-phi)./(1+betta*c_ion));
        M = M_trans+Mconv+M_bc_c+Msc;
        RHS = RHS_trans+RHS_bc_c;
        c_ion = solvePDE(m,M, RHS);
    end
    c_ion_old = c_ion;
    
    % solve the ions diffusion in the water films
    for i = 1:Nx
        BC{i}.left.a = 0; BC{i}.left.b=1; BC{i}.left.c=c_ion.value(i+1); % left boundary
        BC{i}.right.a = 0; BC{i}.right.b=1; BC{i}.right.c=c_ion.value(i+1); % right boundary
        [Mbc, RHSbc] = boundaryCondition(BC{i});
        for j=1:3
            [M_trans, RHS_trans] = transientTerm(c_old{i}, dt, 1.0+a*k_lang./(1+betta*c{i}));
            M = M_trans-Mdiff+Mbc;
            RHS = RHS_trans+RHSbc;
            c{i} = solvePDE(m_diff,M, RHS);
        end
        
        ind_80 = find(c{i}.value(2:end-1)<0.8*c_inj_ion,1, 'first');
        if isempty(ind_80)
            SF.value(i+1) = 1.0;
        else
            SF.value(i+1) = 2*(ind_80-1)/(N-2);
        end
        c_old{i} = c{i};
    end
    
    SF_face = arithmeticMean(SF);
    
    dsw=max(abs(sw_new(:)-sw_old.value(:))./sw_new(:));
    t=t+dt;
    fprintf(1,'\b\b\b\b\b\b\b\b\b\b\b\b\b\bProgress: %d %%',floor(t/t_end*100));
    dt=min([dt*(dsw_alwd/dsw), 1.05*dt, t_end-t]);
    p_old = p;
    sw_old = sw;
    rec_fact=[rec_fact (oil_init-domainInt(1-sw))/oil_init];
    t_day=[t_day t];
    figure(1);visualizeCells(sw); drawnow;
    figure(2); visualizeCells(c_ion);
    
%     figure(2);plot(t_day/3600/24, rec_fact)
%     xlabel('time [day]');
%     ylabel('recovery factor');
%     title([num2str(t/3600/24) ' day']); drawnow;
end
fprintf('\n')
