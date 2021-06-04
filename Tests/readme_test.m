% written by Ali A. Eftekhari
clc
L = 50;  % domain length
Nx = 20; % number of cells
m = createMesh3D(Nx, Nx, Nx, L, L, L);
BC = createBC(m); % all Neumann boundary condition structure
BC.left.a(:) = 0; BC.left.b(:)=1; BC.left.c(:)=1; % Dirichlet for the left boundary
BC.right.a(:) = 0; BC.right.b(:)=1; BC.right.c(:)=0; % right boundary
D_val = 1; % value of the diffusion coefficient
D = createCellVariable(m, D_val); % assign the diffusion coefficient to the cells
D_face = harmonicMean(D); % calculate harmonic average of the diffusion coef on the cell faces
Mdiff = diffusionTerm(D_face); % matrix of coefficients for the diffusion term
[Mbc, RHSbc] = boundaryCondition(BC); % matix of coefficients and RHS vector for the BC
M = Mdiff + Mbc; % matrix of cefficients for the PDE
c = solvePDE(m,M, RHSbc); % send M and RHS to the solver
visualizeCells(c); % visualize the results