function SteadyLidDrivenCavityExample( varargin )
%SteadyLidDrivenCavityExample() solves the steady Lid-driven cavity problem
%   In this 2D example, the momentum equation is solved for velocity and 
%   pressure using a collocated grid with rhie chow interpolation of the 
%   velocies to the cellfaces and the SIMPLE algorithm. The results are
%   then compared to a reference solution given by:
%     Ghia, U. K. N. G., Kirti N. Ghia, and C. T. Shin. 
%     "High-Re solutions for incompressible flow using the Navier-Stokes 
%     equations and a multigrid method."
%     Journal of computational physics 48.3 (1982): 387-411.
%
%   See also MomentumEq, RhieChow, GhiaSolution.
%
%   Date:   18.10.2015
%   Author: Kai Schueller
%
%   Known problems: - works only for equidistant dx dy
%                   - it works for both Rhie Chow and arithmetic mean
%                     interpolation? Maybe checkerboarding effect is not
%                     that big for the considered problem
%                   - it seems that there is no difference when using
%                     arithmetic mean or central difference for Rhie Chow,
%                     but in literature, Rhie Chow needs central difference
%                     to work properly

close all; % closes all figures

addpath('Functions'); % add Functions folder to the search path

% get the reference solution values
testSol=GhiaSolution;

%% User Input
USE_CHOW_INTERP=1;             % flag to choose the rhie-chow interpolation
p_relax=1;                     % pressure relaxation
velo_relax=0.9;                % velocity relaxation
u_init=0;                      % initial guess for x-velocity
v_init=0;                      % initial guess for y-velocity
p_init=0;                      % initial guess for pressure
Re=100;                        % Reynolds number % 100, 400 and 1000 is possible
lidVelo=1;                     % velocity of the lid [m/s]
rho=1000;                      % density [kg/m^3]
cavityLength=0.1;              % length of one edge of the quare cavity [m]
mu=rho*cavityLength*lidVelo/Re;% dynamic viscosity
nIter=250;                     % number of iterations

% First create the mesh
m=createMesh2D(51,51,cavityLength,cavityLength);
%visualizeMesh2D(m); % Uncomment if mesh should be displayed

%% Boundary Conditions
% Since the momentum equation is solved for both, u (velocity in x) and v
% (velocity in y), we need two separate sets of boundaries (u and v)
U_BC = createBC(m);
V_BC = createBC(m);

% Assign boundary conditions u-momentum equation
U_BC.top.a(:) = 0; U_BC.top.b(:)=1; U_BC.top.c(:)=lidVelo;    % lid velocity in x-direction at the top of the cavity
U_BC.left.a(:) = 0; U_BC.left.b(:)=1; U_BC.left.c(:)=0;       % zero x-velocity at left boundary 
U_BC.right.a(:) = 0; U_BC.right.b(:)=1; U_BC.right.c(:)=0;    % zero x-velocity right boundary 
U_BC.bottom.a(:) = 0; U_BC.bottom.b(:)=1; U_BC.bottom.c(:)=0; % no-slip at bottom boundary 

% Assign boundary conditions for v-momentum equation
V_BC.top.a(:) = 0; V_BC.top.b(:)=1; V_BC.top.c(:)=0;          % zero y-velocity at top boundary
V_BC.bottom.a(:) = 0; V_BC.bottom.b(:)=1; V_BC.bottom.c(:)=0; % zero y-velocity at bottom boundary
V_BC.left.a(:) = 0; V_BC.left.b(:)=1; V_BC.left.c(:)=0;       % no-slip at left boundary 
V_BC.right.a(:) = 0; V_BC.right.b(:)=1; V_BC.right.c(:)=0;    % no-slip at right boundary 

% Assign boundary conditions for the pressure correction equation
P_BC = createBC(m); % this automatically creates all Neumann, which is needed here
% we need to set a refence pressure, because otherwise matlab throws the
% following warning for pressure correction: Matrix is close to singular or
% badly scaled. Results may be inaccurate.
P_BC.bottom.a(end/2+0.5) = 0; P_BC.bottom.b(end/2+0.5)=1; P_BC.bottom.c(end/2+0.5)=0;

%% build required cell and face variables using initital guess for velocities and pressure
mu_faceVar=createFaceVariable(m,mu);

rho_faceVar=createFaceVariable(m,rho);

U_cellVar=createCellVariable(m,u_init,U_BC);
V_cellVar=createCellVariable(m,v_init,V_BC);

faceVelocity=createFaceVariable(m,[u_init v_init]);

p_cellVar=createCellVariable(m,p_init,P_BC);

pGradxCellVar=createCellVariable(m,0);
pGradyCellVar=pGradxCellVar;

%% main loop
for i=1:nIter
    %% Momentum equation
    %pGradMat = gradientCellTerm(p_cellVar);
    
    % for Rhie Chow, we need central difference
    pGradMat = CDGradientCellTerm(p_cellVar);
    
    U_cellVar_old=U_cellVar;
    V_cellVar_old=V_cellVar;

    % solve the momentum equations
    [U_cellVar,V_cellVar,U_aP_sumMat,V_aP_sumMat]=MomentumEq(velo_relax,m,pGradMat,mu_faceVar,faceVelocity,U_cellVar_old,V_cellVar_old,U_BC,V_BC,rho,'SIMPLE');
        
    %% Pressure Correction

    % interpolate velocity at faces using rhie chow interpolation
    pGradxCellVar.value=pGradxCellVar.value*0;
    pGradyCellVar.value=pGradyCellVar.value*0;
    pGradxCellVar.value(2:end-1,2:end-1)=pGradMat.xvalue;
    pGradyCellVar.value(2:end-1,2:end-1)=pGradMat.yvalue;
    
    if USE_CHOW_INTERP
        faceVelocity=RhieChow(U_cellVar,V_cellVar,U_aP_sumMat,V_aP_sumMat,p_cellVar,pGradxCellVar,pGradyCellVar,U_BC,V_BC);
    else
        helpx=arithmeticMean(U_cellVar);
        helpy=arithmeticMean(V_cellVar);
        faceVelocity.xvalue=helpx.xvalue;
        faceVelocity.yvalue=helpy.yvalue;
    end

    % ---------------------------------------------------------------------
    % calculate divergence of the velocity
    [RHS, ~, ~, ~] = divergenceTerm(faceVelocity);
    % ---------------------------------------------------------------------

    % ---------------------------------------------------------------------
    % calculate the diffusion term
    U_aP_Face=arithmeticMean(U_aP_sumMat);
    V_aP_Face=arithmeticMean(V_aP_sumMat);
    diffusion_coeff_helpVar=U_aP_Face;
    diffusion_coeff_helpVar.yvalue=V_aP_Face.yvalue;

    diffusion_coeff.domain=m;
    diffusion_coeff.xvalue=velo_relax*rho_faceVar.domain.cellsize.y(1)./diffusion_coeff_helpVar.xvalue;
    diffusion_coeff.yvalue=velo_relax*rho_faceVar.domain.cellsize.x(1)./diffusion_coeff_helpVar.yvalue;

    test_diff=diffusionTerm(diffusion_coeff);
    % ---------------------------------------------------------------------
    
    % ---------------------------------------------------------------------
    % pressure boundary condition
    [p_Mbc, p_RHSbc] = boundaryCondition2D(P_BC);
    % ---------------------------------------------------------------------

    % ---------------------------------------------------------------------
    % solve for the new pressure correction
    p_new=solvePDE(m,test_diff+p_Mbc, RHS+p_RHSbc);
    % ---------------------------------------------------------------------
    
    p_cellVarOld.value=p_cellVar.value;
    p_cellVar.value=p_cellVar.value+p_relax*p_new.value;
    p_max_error=max(max(abs(p_cellVarOld.value-p_cellVar.value)));

    % ---------------------------------------------------------------------
    % update the velocity
    pGradMat = gradientCellTerm(p_new);
    pGradxCellVar.value=pGradxCellVar.value*0;
    pGradyCellVar.value=pGradyCellVar.value*0;
    pGradxCellVar.value(2:end-1,2:end-1)=pGradMat.xvalue;
    pGradyCellVar.value(2:end-1,2:end-1)=pGradMat.yvalue;
    
    UU=U_cellVar.value-...
        pGradxCellVar.domain.cellsize.y(1)*pGradxCellVar.value./U_aP_sumMat.value;
    VV=V_cellVar.value-...
        pGradxCellVar.domain.cellsize.x(1)*pGradyCellVar.value./V_aP_sumMat.value;

    U_cellVar.value=velo_relax*UU+(1-velo_relax)*U_cellVar.value;
    V_cellVar.value=velo_relax*VV+(1-velo_relax)*V_cellVar.value;
    % ---------------------------------------------------------------------

    % plotting at every iteration would be very slow, so we just plot at
    % every 50th iteration
    if ~mod(i,50)
        contour(m.cellcenters.x,m.cellcenters.y,hypot(U_cellVar.value(2:end-1,2:end-1)',V_cellVar.value(2:end-1,2:end-1)'),'LevelList',0:0.1:1,'ShowText','on')
        hold on
        %quiver(m.cellcenters.x,m.cellcenters.y,U_cellVar.value(2:end-1,2:end-1)',V_cellVar.value(2:end-1,2:end-1)');
        axis equal
        xlabel('x [m]');
        ylabel('y [m]');
        hold off
        ylim([0 0.1]);
        xlim([0 0.1]);

        drawnow
    end
    
    %% Calculate Root Mean Squared Error
    % First, we need to inter/extrapolate the datapoints to the datapoints
    % given by the reference solution
    GhiaValues.u=interp1(m.cellcenters.y,U_cellVar.value(end/2+0.5,2:end-1),testSol.y*0.1,'linear','extrap');
    GhiaValues.v=interp1(m.cellcenters.x,V_cellVar.value(end/2+0.5,2:end-1),testSol.x*0.1,'linear','extrap');
    
    % Now, we calculate the root mean squared error for u and v velocities
    RMSEu=sqrt(sum((testSol.u{testSol.Re==Re}(:)-GhiaValues.u(:)).^2)/numel(testSol.u{testSol.Re==Re}));
    RMSEv=sqrt(sum((testSol.v{testSol.Re==Re}(:)-GhiaValues.v(:)).^2)/numel(testSol.v{testSol.Re==Re}));
    
    U_max_error=max(max(abs(U_cellVar_old.value-U_cellVar.value)));
    V_max_error=max(max(abs(V_cellVar_old.value-V_cellVar.value)));
    
    fprintf('i: %d\t p_max_error: %e\t (u, v)_max_error: (%e, %e)\t RMSE(u, v): (%e, %e)\n',i,p_max_error,U_max_error,V_max_error,RMSEu,RMSEv);
end

%% Plot the final solution versus the reference solution
figure
hold on;
h1=plot(testSol.u{testSol.Re==Re},testSol.y*0.1,'o');
h2=plot(U_cellVar.value(end/2+0.5,2:end-1),m.cellcenters.y);
legend([h1 h2],'Ghia et al. [1]','Numerical solution','location','southeast');
ylabel('y [m]')
xlabel('Velocity in x-direction [m/s]');
hold off;

figure;
hold on;
h1=plot(testSol.x*0.1,testSol.v{testSol.Re==Re},'o');
h2=plot(m.cellcenters.x,V_cellVar.value(2:end-1,end/2+0.5));
legend([h1 h2],'Ghia et al. [1]','Numerical solution');
xlabel('x [m]')
ylabel('Velocity in y-direction [m/s]');
hold off;
end

function cellGrad=CDGradientCellTerm(phi)
	Nx = phi.domain.dims(1);
	Ny = phi.domain.dims(2);
	DX = repmat(phi.domain.cellsize.x(2:end-1), 1, Ny);
	DY = repmat(phi.domain.cellsize.y(2:end-1)', Nx, 1);
	xvalue = (phi.value(3:Nx+2,:)-phi.value(1:Nx,:))./(2*DX(2));
	yvalue = (phi.value(:,3:Ny+2)-phi.value(:,1:Ny))./(2*DY(2));
	zvalue=[];
	cellGrad=FaceVariable(phi.domain, xvalue(:,2:end-1), yvalue(2:end-1,:), zvalue);
end