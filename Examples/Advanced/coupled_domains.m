% Coupling two domains: this is my first test and I'm not sure if it really
% works
% heat transfer between two 1D domains; it can be solved without copling in
% only one domain, which helps to compare the results
% tried to couple two domainsl does not work as far as I can tell :-(
% Written by Ali A. Eftekhari
% Last checked: June 2021
Nx = 10;
Lx = 0.2;
k1 = 0.1; % J/(m.s.K)
k2 = 0.5; % J/(m.s.K)
T_left = 300; % K
T_right = 350; % K

%% solution on a single domain
m = createMesh1D(2*Nx, 2*Lx);
k = [k1*ones(Nx, 1); k2*ones(Nx,1)];
k_cell = createCellVariable(m, k);
k_face = harmonicMean(k_cell);
BC = createBC(m);
BC.left.a(:) = 0.0; BC.left.b(:) = 1.0; BC.left.c(:) = T_left;
BC.right.a(:) = 0.0; BC.right.b(:) = 1.0; BC.right.c(:) = T_right;
[M_bc, RHS_bc] = boundaryCondition(BC);
M_diff = diffusionTerm(k_face);
T = solvePDE(m, M_bc-M_diff, RHS_bc);
visualizeCells(T);

%% solution on two coupled domains (does not work)
m1 = createMesh1D(Nx, Lx);
m2 = createMesh1D(Nx, Lx);
m2.cellcenters.x = m2.cellcenters.x + Lx;
m2.facecenters.x = m2.facecenters.x + Lx;

k1_cell = createCellVariable(m1, k1);
k2_cell = createCellVariable(m2, k2);
k1_face = harmonicMean(k1_cell);
k2_face = harmonicMean(k2_cell);
M_diff_1 = diffusionTerm(k1_face);
M_diff_2 = diffusionTerm(k2_face);

% now the interesting part: boundary conditions
% we have a constant flux from the right boundary of the first domain to
% the left boundary of the second domain
BC1_nc = createBC(m1);
BC2_c = createBC(m2);
BC1_c = createBC(m1);
BC2_nc = createBC(m2);
% not-coupled boundaries:
BC1_nc.left.a(:) = 0.0; BC1_nc.left.b(:) = 1.0; BC1_nc.left.c(:) = T_left;
BC2_nc.right.a(:) = 0.0; BC2_nc.right.b(:) = 1.0; BC2_nc.right.c(:) = T_right; 
% clean the rest (I will add a convenience function for this, later):
BC1_nc.right.a(:) = 0.0; BC1_nc.right.b(:) = 0.0; BC1_nc.right.c(:) = 0.0;
BC2_nc.left.a(:) = 0.0; BC2_nc.left.b(:) = 0.0; BC2_nc.left.c(:) = 0.0;

% coupled boundaries:
BC1_c.right.a(:) = k1_face.xvalue(end); BC1_c.right.b(:) = 0.0; BC1_c.right.c(:) = 0.0;
BC2_c.left.a(:) = k2_face.xvalue(1); BC2_c.left.b(:) = 0.0; BC2_c.left.c(:) = 0.0;
% clean the rest:
BC1_c.left.a(:) = 0.0; BC1_c.left.b(:) = 0.0; BC1_c.left.c(:) = 0.0;
BC2_c.right.a(:) = 0.0; BC2_c.right.b(:) = 0.0; BC2_c.right.c(:) = 0.0;

[M_bc_1_c, RHS_1_c] = boundaryCondition(BC1_c);
[M_bc_2_c, RHS_2_c] = boundaryCondition(BC2_c);
[M_bc_1_nc, RHS_1_nc] = boundaryCondition(BC1_nc);
[M_bc_2_nc, RHS_2_nc] = boundaryCondition(BC2_nc);

M = [-M_diff_1+M_bc_1_c+M_bc_1_nc -M_bc_2_c
    M_bc_1_c -M_diff_2+M_bc_2_nc-M_bc_2_c];
RHS = [RHS_1_c-RHS_2_c+RHS_1_nc; RHS_1_c-RHS_2_c+RHS_2_nc];

T_new = M\RHS;
hold on
plot(T_new)
