% another simple convection example using method of lines, still a bit
% messy, wait for a blog post
clc; clear;
% define a 2D mesh
H = 1;
W = 1;
Nx = 100;
Ny = 100;
mesh1 = buildMesh2D(Nx, Ny, W, H);
% velocity field
u = 0.001*ones(Nx,Ny);
uf = arithmeticMean(mesh1, u);
% uf.yvalue(:,:)=0;
% diffusion field
D = 1e-2*createCellVariable(mesh1, 1e-2);
Df = arithmeticMean(mesh1, D);
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
phi.Old = createCellVariable(mesh1, 0, BC);
phi.Old(5:10, 5:10) = 1;
% define the convection term
Mconv = convectionUpwindTerm(mesh1, uf);
Mdif = diffusionTerm(mesh1, Df);
% define the BC term
[Mbc, RHSbc] = boundaryCondition(mesh1, BC);
% solver
final_t = 500;
% define the transient term
[M, RHS] = combineBC2D(mesh1, BC, Mconv, ...
    zeros((Nx+2)*(Ny+2),1));
% eq: dcdt = D d2c/dx2
dcdt = @(t,c)(-M*c-RHS);
[t_temp, c_temp] = ode45(dcdt, [0 final_t], phi.Old(:));
phi.value = reshape(c_temp(end,:), Nx, Ny);
figure(1);pcolor(phi.value');colorbar;
