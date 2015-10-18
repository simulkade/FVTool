function CouetteFlow( varargin )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
close all; % closes all figures
addpath('../Functions');

h=4; % height of the channel
U=10; % velocity of the moving surface
dPdx=0; % pressure gradient in x direction
mu=0.001; % dynamic viscosity

% function for the analytical solution
u_analFunc=@(y) y./h.*(U-h.^2./(2*mu).*dPdx.*(1-y/h));

% create mesh
m=createMesh2D(8,8,2*h,h);

%% Boundary Conditions
U_BC = createBC(m); % all Neumann boundary condition structure 

% Assign boundary conditions
U_BC.left.a(:) = 1; U_BC.left.b(:)=0; U_BC.left.c(:)=0; % Dirichlet for the left boundary 
U_BC.right.a(:) = 1; U_BC.right.b(:)=0; U_BC.right.c(:)=0; % Dirichlet for the left boundary 
U_BC.top.a(:) = 0; U_BC.top.b(:)=1; U_BC.top.c(:)=U; % Dirichlet for the left boundary 
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

mu_faceVar=createFaceVariable(m,0.001);

rho=1000;

pGradMat = gradientCellTerm(p_cellVar);

faceVelocity=createFaceVariable(m,[0 0]);
U_cellVar=createCellVariable(m,0,U_BC);
V_cellVar=createCellVariable(m,0,V_BC);

% solve momentum equation Momentum equation
[U_cellVar,~,~,~]=MomentumEq(1,m,pGradMat,mu_faceVar,faceVelocity,U_cellVar,V_cellVar,U_BC,V_BC,rho,'SIMPLE');

y_anal=linspace(0,h,100);
u_anal=u_analFunc(y_anal);

figure;
plot(u_anal,y_anal); % plot analytical solution
hold on;
plot(U_cellVar.value(3,2:end-1),U_cellVar.domain.cellcenters.y,'o');
hold off;
figure;
visualizeCells(U_cellVar);

end

