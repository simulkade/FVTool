% another simple convection example using method of lines, still a bit
% messy, wait for a blog post
% Using Matlab ODE solvers to solve a convective flux that 
% moves a plume in diagonal direction
% see how terribly diffusive the solution is!
% Written by Ali A. Eftekhari
% Last checked: June 2021
clc
% define a 2D mesh
H = 1;
W = 1;
Nx = 100;
Ny = 100;
mesh1 = createMesh2D(Nx, Ny, W, H);
% velocity field
u = createCellVariable(mesh1, 0.001*ones(Nx,Ny));
uf = arithmeticMean(u);
% uf.yvalue(:,:)=0;
% diffusion field
D = 1e-2*createCellVariable(mesh1, 1e-2);
Df = arithmeticMean(D);
% transient term coefficient
alfa = createCellVariable(mesh1, 1);
% define the boundaries
% dirichlet on all the boundaries
BC = createBC(mesh1); % all Neumann
BC.bottom.a(:) = 0; BC.bottom.b(:) = 1; BC.bottom.c(:) = 0;
BC.left.a(:) = 0; BC.left.b(:) = 1; BC.left.c(:) = 0;
% BC.right.a(:) = 0; BC.right.b(:) = 1; BC.right.c(:) = 0;
% BC.top.a(:) = 0; BC.top.b(:) = 1; BC.top.c(:) = 0;
% Initial values
phi_old = createCellVariable(mesh1, 0, BC);
phi_old.value(5:10, 5:10) = 1;
% define the convection term
Mconv = convectionUpwindTerm(uf);
Mdif = diffusionTerm(Df);
% define the BC term
[Mbc, RHSbc] = boundaryCondition(BC);
% solver
final_t = 500;
% define the transient term
[M, RHS] = combineBC2D(BC, Mconv, ...
    zeros((Nx+2)*(Ny+2),1));
% eq: dcdt = D d2c/dx2
dcdt = @(t,c)(-M*c-RHS);
[t_temp, c_temp] = ode45(dcdt, [0 final_t], internalCells(phi_old));
for i =1:length(t_temp)
    phi.value = reshape(c_temp(i,:), Nx, Ny);
    figure(1);pcolor(phi.value');title(["t = ", num2str(t_temp(i))]); colorbar;drawnow;
end