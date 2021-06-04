% Written by Ali A. Eftekhari
% Last checked: June 2021
clc % clean the command prompt
L = 0.1; % 10 cm length
W = 0.01; % 1 cm thickness
Nx = 50; % number of cell in x direction
Ny = 20; % number of cells in the y direction
m = createMesh2D(Nx, Ny, L, W); % creates a 2D Cartesian grid
% m = createMesh3D(Nx, Ny,Ny, L, W,W); % creates a 3D Cartesian grid, activate this line if you are curious
T_inf = 25+273.15; % [K] ambient temperature
T_base = 100+273.15; % [K] temperature at the base of the fin
k_val = 237; % W/(m.K) thermal conductivity
h_val = 10; % W/(m^2.K) heat transfer coefficient
k = createCellVariable(m, k_val); % assign thermal cond. value to all cells
k_face = geometricMean(k); % geometric average of the thermal conductivity values on the cell faces
BC = createBC(m); % creates a BC structure for the domain m; all Neumann boundaries
BC.left.a(:)=0; BC.left.b(:)=1; BC.left.c(:)=T_base; % convert the left boundary to constant temperature
BC.right.a(:)=k_val/h_val; BC.right.b(:)=1; BC.right.c(:)=T_inf; % right boundary to Robin
BC.top.a(:)=k_val/h_val; BC.top.b(:)=1; BC.top.c(:)=T_inf; % top boundary to Robin
BC.bottom.a(:)=-k_val/h_val; BC.bottom.b(:)=1; BC.bottom.c(:)=T_inf; % bottom boundary to Robin (don't forget the normal vector sign)
M_cond = diffusionTerm(k_face); % matrix of coefficients for the heat diffusion
[M_bc, RHS_bc] = boundaryCondition(BC); % matrix of coefficients and RHS vector for the boundary conditions
T = solvePDE(m, M_cond+M_bc, RHS_bc); % solve the linear system of discretized PDE
visualizeCells(T); % visualize the results
