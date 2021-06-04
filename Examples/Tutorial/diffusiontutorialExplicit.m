% a tutorial adapted from the fipy diffusion 1D example
% see: http://www.ctcms.nist.gov/fipy/examples/diffusion/index.html

% Written by Ali A. Eftekhari
% Last checked: June 2021
clc

%% define the domain
L = 5;  % domain length
Nx = 100; % number of cells
meshstruct = createMesh1D(Nx, L);
BC = createBC(meshstruct); % all Neumann boundary condition structure
BC.left.a(:) = 0; BC.left.b(:)=1; BC.left.c(:)=1; % left boundary
BC.right.a(:) = 0; BC.right.b(:)=1; BC.right.c(:)=0; % right boundary
x = meshstruct.cellcenters.x;
%% define the transfer coeffs
D_val = 1;
alfa = createCellVariable(meshstruct, 1);
Dave = createFaceVariable(meshstruct, D_val);
%% define initial values
c_old = createCellVariable(meshstruct, 0, BC); % initial values
c = c_old;
%% loop
dt = 0.001; % time step
final_t = 0.5;
for t=dt:dt:final_t
    % step 1: calculate divergence term
    RHS = divergenceTerm(Dave.*gradientTerm(c_old));
    % step 2: calculate the new value for internal cells
    c_old_int=internalCells(c_old);
    c_val_int = dt*excludeGhostRHS(meshstruct, RHS+constantSourceTerm(c_old))+c_old_int(:);
    c_analytical = 1-erf(x/(2*sqrt(D_val*t)));
    figure(1);plot(x, c_val_int, x, c_analytical, 'r--');drawnow;
    % Step 3: calculate the ghost cell values
    c.value=cellBoundary(c_val_int, BC);
    c_old = c;
    disp(t)
end
