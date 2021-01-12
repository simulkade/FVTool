%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Version 1, two dimension, Panga Model(2006)           %%%%
%%%% Developed by Behzad Hosseinzadeh                      %%%%
%%%% Two-scale continuum model, more information:          %%%%
%%%% https://doi.org/10.1002/aic.10574                     %%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; clear;
%% define the domain
nx      = 100;
L       = 15.24;                % Core Length [cm]
H       = 3.81;                 % Core width [cm]
h       = 1.0;                  % Dimensionless Core width [cm]
alpha0  = H / L;                % Aspect Ratio
ny      = alpha0 * nx;          % Calculate ny according to nx
r0      = 0.001;                % Initial mean pore size [cm]
a0      = 50.000;               % Initial interfacial area per unit volume [cm^-1]
K0      = 0.100;                % Initial average permeability [D]
ks      = 0.2;                  % Surface dissolution reaction-rate constant [cm/s]
Dm0     = 3.60e-5;              % Acid diffusivity [cm^2/s]
MWofA   = 36.46;                % molecular weight of Acid (gr / mol)
pHin    = 0.01; 
rhoSolid= 2.71;                 % Rock density [g/cm^3]
nho     = 0.01;                 % Fluid  kinematic  viscosity ...
                                % Kinematic viscosity, is the ratio of the ...
                                % dynamic viscosity to the density. The dimensions of ...
                                % kinematic viscosity are L2 [cm^2/s].
pExit   = 68;                   % Exit pressure [atm]
mio0    = 1;                    % Acid viscosity [cp]
ShAsym  = 3.66;                 % Asymptotic Sherwood number
eps0    = 0.20;                 % Average porosity
delEps  = 0.15;                 % Heterogeneity magnitude
lHT     = 1;                    % Heterogeneity length scale [mm]
beta    = 2.00;                 % Pore-broadening parameter
gama    = 1;                    % Pore-connectivity parameter
alpha0s = 0.5;                  % Constant in dispersion correlations
landaX  = 0.5;                  % Constant in axial dispersion correlation
landaT  = 0.1;                  % Constant in transverse dispersion correlation
mp      = 1.0;                  % Constant in mass-transfer correlation 
u0      = 0.01;                 % The injection velocities of the acid [0.5 MHCl] are varied between 0.14 [cm/s] ...
                                % and 1.4e-4 [cm/s], where 0.14 [cm/s] corresponds to the uniform dissolution ...
                                % regime and 1.4e-4 (cm/s) corresponds to the face dissolution regime
%%%%%%%% dt correction %%%%%%%%%
dt      = 0.003 * L / u0;       % Time [s]
%%%%%%%% dt correction %%%%%%%%%
tTotal  = L / u0;               % Total Time
BTtime  = 0;                    % BreakThrough time
%% Calculate initial Variables
Cf0     = 0.05;                 % Inlet  Acid  Concentration
rhoAcid = 1.00;                 % Acid density [g/cm^3]
rhoWater= 1;                    % Water density [g/cm^3]
Cf0     = Cf0 * (Cf0 * rhoAcid + (1 - Cf0) * rhoWater);
alpha   = 1.370;                % Acid  Dissolving  Power 
Sch     = nho / Dm0;            % Schmidt Number 
%% Functions 
kCal        = @(Eps)(K0 .* power((Eps ./ eps0), gama) .* power((Eps .* (1-eps0))./...
               (eps0 .* (1-Eps)), 2*beta));
rCal        = @(Eps,K)(r0 .* power((K .* eps0)./(K0 .* Eps),1/2));
AvCal       = @(Eps,rp)(a0 .* (Eps .* r0) ./ (eps0 .* rp));
ReCal       = @(u,rp)(2 .* rp .* u ./ nho);
ShCal       = @(Re)(ShAsym + 0.7 .* power(Re,1/2) .* power(Sch,1/3) ./ power(mp,1/2));
kcCal       = @(Sh,Dm,rp)(Sh .* Dm ./ 2 ./ rp);
dpCoreCal   = @(p)(sum(p.value(1,2:ny+1)+p.value(2,2:ny+1))/ny/2 ...
            - sum(p.value(nx+1,2:ny+1)+p.value(nx+2,2:ny+1))/ny/2);
DifCalX     = @(Eps,u_center,rp,Dm)(alpha0s .* Dm .* Eps + 2 .* landaX .* u_center .* rp);
DifCalT     = @(Eps,u_center,rp,Dm)(alpha0s .* Dm .* Eps + 2 .* landaT .* u_center .* rp);
% This function returns pH without logarithm, -log10() need to be run then
pHCal       = @(cf)((cf .* rhoAcid) .* 1000 ./ (MWofA));
%% Initialization 
% create mesh and assign values to the domain
m   = createMesh2D(nx, ny, L, H);
% create mesh and assign values to the domain
% Initialize Porosity
Eps(:,:) = eps0 + delEps + (-delEps-delEps).*rand(nx,ny);
% Initialize variable
Eps     = createCellVariable(m, Eps);
K       = kCal(Eps);
rp      = rCal(Eps,K);
av      = AvCal(Eps,rp);
Rep     = ReCal(u0,rp);
Sh      = ShCal(Rep);
kc      = kcCal(Sh,Dm0,rp);
cf      = createCellVariable(m, zeros(nx,ny));
p       = createCellVariable(m, zeros(nx,ny));
pH      = pHCal(cf);
mio     = createCellVariable(m, mio0);
Dmf     = createCellVariable(m, Dm0);
D_face  = createFaceVariable(m, zeros(nx,ny));

itr     = 1;
dpBT    = 0.01;

% Define the boundaries
% a??.n+b?=c. 
BCp = createBC(m); % Neumann BC for pressure
BCc = createBC(m); % Neumann BC for concentration
% change the right boandary to constant pressure (Dirichlet)
BCp.right.a(:) = 0.0;
BCp.right.b(:) = 1.0;
BCp.right.c(:) = pExit;
% left boundary
BCp.left.b(:)  = 0.0;
BCp.left.c(:)  = u0;
% change the left boundary (Robin boundary conditions)
BCc.left.b(:) = u0;
BCc.left.c(:) = u0*Cf0;
while 1
    tic
    % **Pressure**
    perm_face       = harmonicMean(K ./ mio);
    % Set left boundary condition
    BCp.left.a(:)   = - perm_face.xvalue(1,:);
    Mdiffp          = - diffusionTerm(perm_face);
    RHSreaction     = - constantSourceTerm(cf.*av.*alpha.*kc.*ks./(ks+kc)./rhoSolid);
    [Mbcp, RHSb]    = boundaryCondition(BCp);
    Mp              = Mdiffp + Mbcp;
    RHSp            = RHSb + RHSreaction;
    p               = solvePDE(m, Mp, RHSp);
    
    % Determine breakthrough condition
    if itr == 1
        dpBT = dpCoreCal(p) / 100;
    end
    u               = - perm_face .* gradientTerm(p);
    u_center        = power(power((u.xvalue(1:nx,:) + u.xvalue(2:nx+1,:)) / 2,2) + ...
                      power((u.yvalue(:,1:ny) + u.yvalue(:,2:ny+1)) / 2,2),1/2);
    u_center        = createCellVariable(m, u_center);
    
    % **Acid Concentration**
    D_face_temp     = harmonicMean(DifCalX(Eps,u_center,rp,Dmf));
    D_face.xvalue   = D_face_temp.xvalue;
    D_face_temp     = harmonicMean(DifCalT(Eps,u_center,rp,Dmf));
    D_face.yvalue   = D_face_temp.yvalue;
    % Set left boundary condition
    BCc.left.a(:)   = - D_face.xvalue(1,:);
    BCc.left.c(:)   = u0*Cf0;
    Mdiffc          = - diffusionTerm(D_face);
    Mconvc          = convectionUpwindTerm(u);
    Mlinc           = linearSourceTerm(kc.*ks.*av./(kc+ks));
%     Mcontc          = - continuityVelocityTerm(u);
    [Mbcc, RHSbcc]       = boundaryCondition(BCc);
    [Mtransc, RHStransc] = transientTerm(cf, dt, Eps);
    RHSc            = RHSbcc+RHStransc;
    Mc              = Mdiffc+Mlinc+Mconvc+Mbcc+Mtransc;
    cf              = solvePDE(m, Mc, RHSc);
    cf.value(cf.value>1) = Cf0;
    
    % Update local porosity 
    Rep             = ReCal(u_center,rp);
    Sh              = ShCal(Rep);
    kc              = kcCal(Sh,Dmf,rp);
    Eps_old         = Eps;
    Eps             = (cf.*av.*alpha.*kc.*ks./(ks+kc)./rhoSolid).*dt+Eps;
    Eps.value(1,:)  = 0.99;
    
    % Update pore scale calculation
    K               = kCal(Eps);
    rp              = rCal(Eps,K);
    av              = AvCal(Eps,rp);
    Rep             = ReCal(u_center,rp);
    Sh              = ShCal(Rep);
    kc              = kcCal(Sh,Dmf,rp);
    pH              = pHCal(cf);
    pH.value(:)     = - log10(pH.value(:));
    pH.value(pH.value > 6.99) = 6.99;
    pH.value(pH.value < pHin) = pHin;
    
    if dpCoreCal(p) < dpBT
        break;
    end
    
    if rem(itr,10) == 0
        clc;
        figure(1);visualizeCells(Eps);drawnow;
    end
    
%     if rem(itr,500) == 0 || itr == 1
%         save(strcat('Eps',int2str(itr)),'Eps');
%     end
    
    BTtime          = BTtime + dt;
    itr = itr + 1;
    toc
end
PVBT = u0 * BTtime / eps0 / L
figure(1);visualizeCells(Eps);drawnow;
% save(strcat('Eps',int2str(itr)),'Eps');