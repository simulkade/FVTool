function PoiseuilleFlow( varargin )
%POISOILLEFLOW Summary of this function goes here
%   Detailed explanation goes here
close all; % closes all figures
addpath('../Functions');

h=4; % height of the channel
dPdx=-0.01; % pressure gradient in x direction
mu=0.001; % dynamic viscosity

% function for the analytical solution
u_analFunc=@(y) -0.5*h^2/mu*dPdx*(1-(y/h).^2);

% create mesh
m=createMesh2D(20,20,2*h,2*h);

%% Boundary Conditions
U_BC = createBC(m); % all Neumann boundary condition structure 

% Assign boundary conditions
U_BC.left.a(:) = 1; U_BC.left.b(:)=0; U_BC.left.c(:)=0; % Dirichlet for the left boundary 
U_BC.right.a(:) = 1; U_BC.right.b(:)=0; U_BC.right.c(:)=0; % Dirichlet for the left boundary 
U_BC.top.a(:) = 0; U_BC.top.b(:)=1; U_BC.top.c(:)=0; % Dirichlet for the left boundary 
U_BC.bottom.a(:) = 0; U_BC.bottom.b(:)=1; U_BC.bottom.c(:)=0; % Dirichlet for the left boundary 

V_BC = createBC(m); % all Neumann boundary condition structure 

% Assign boundary conditions
V_BC.top.a(:) = 0; V_BC.top.b(:)=1; V_BC.top.c(:)=0; % Dirichlet for the left boundary 
V_BC.bottom.a(:) = 0; V_BC.bottom.b(:)=1; V_BC.bottom.c(:)=0; % Dirichlet for the left boundary 
V_BC.left.a(:) = 0; V_BC.left.b(:)=1; V_BC.left.c(:)=0; % Dirichlet for the left boundary 
V_BC.right.a(:) = 0; V_BC.right.b(:)=1; V_BC.right.c(:)=0; % Dirichlet for the left boundary 

BC_pressureCorrection = createBC(m);
BC_pressureCorrection.left.a(:) = 1; BC_pressureCorrection.left.b(:)=0; BC_pressureCorrection.left.c(:)=0;
BC_pressureCorrection.right.a(:) = 1; BC_pressureCorrection.right.b(:)=0; BC_pressureCorrection.right.c(:)=0;
BC_pressureCorrection.bottom.a(:) = 1; BC_pressureCorrection.bottom.b(:)=0; BC_pressureCorrection.bottom.c(:)=0;
BC_pressureCorrection.bottom.a(end/2) = 0; BC_pressureCorrection.bottom.b(end/2)=1; BC_pressureCorrection.bottom.c(end/2)=1;
BC_pressureCorrection.top.a(:) = 1; BC_pressureCorrection.top.b(:)=0; BC_pressureCorrection.top.c(:)=0;

% constant pressure - so there is no gradient
p_cellVar=createCellVariable(m,1,BC_pressureCorrection);

mu_cellVar=createCellVariable(m,0.001);
mu_faceVar=createFaceVariable(m,0.001);

rho=1000;
rho_cellVar=createCellVariable(m,rho);
rho_faceVar=createFaceVariable(m,rho);

U_cellVar=createCellVariable(m,0,U_BC);
V_cellVar=createCellVariable(m,0,V_BC);
faceVelocity=createFaceVariable(m,[0 0]);

pGradMat = gradientCellTerm(p_cellVar);
helpVarx = createCellVariable(m,dPdx);
helpVary = createCellVariable(m,0);
pGradMat.xvalue=helpVarx.value(2:end-1,2:end-1);
pGradMat.yvalue=helpVary.value(2:end-1,2:end-1);

% solve momentum equation Momentum equation
[U_cellVar,~,~,~]=MomentumEq(1,m,pGradMat,mu_faceVar,faceVelocity,U_cellVar,V_cellVar,U_BC,V_BC,rho,'SIMPLE');

y_anal=linspace(0,h,100);
u_anal=u_analFunc(y_anal);

figure;
plot(u_anal,y_anal); % plot analytical solution
hold on;
plot(U_cellVar.value(3,2:end-1),U_cellVar.domain.cellcenters.y-h,'o');
hold off;
figure;
visualizeCells(U_cellVar);


end

