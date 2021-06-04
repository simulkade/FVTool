%% Transient diffusion equation
% This example solves diffusion-adsorption equation
% assume a fluid domain with fixed concentration on both ends. A solute
% diffuses into the domain and gets adsorbed on a solid surface.
% If the concentration in the domain goes above 80% of the boundary 
% concentration, we record it as a change of surface properties and finally plot how the 
% surface properties of the whole domain changes versus time.
% this is part of our paper:
% https://pubs.acs.org/doi/abs/10.1021/acs.energyfuels.0c02493
% Written by Ali A. Eftekhari
% Last checked: June 2021
clc
%% partameters
a = 1e9; % m^2/m^3
k = 1e-4; % Langmuir adsorption coefficient
betta = 1.0e-4; % Langmuir adsorption coefficient 2
%% Define the domain and create a mesh structure
L = 1e-6; % [m]  % domain length
Nx = 50; % number of cells
m = createMesh1D(Nx, L);
%% Create the boundary condition structure
BC = createBC(m); % all Neumann boundary condition structure
BC.left.a = 0; BC.left.b=1; BC.left.c=1; % left boundary
BC.right.a = 0; BC.right.b=1; BC.right.c=1; % right boundary
x = m.cellcenters.x;
%% define the transfer coeffs
D_val = 1e-13; % m^2/s effective diffusivity
D = createCellVariable(m, D_val);
alfa = createCellVariable(m, 1);
%% define initial values
c_init = 0;
c_old = createCellVariable(m, c_init, BC); % initial values
c = c_old; % assign the old value of the cells to the current values
%% loop
Dave = harmonicMean(D);
Mdiff = diffusionTerm(Dave);
[Mbc, RHSbc] = boundaryCondition(BC);
dt = 100; % time step
final_t = 2000*dt;
t_mod_wet = 0:dt:final_t;
mod_area = zeros(size(t_mod_wet));
% hold all
i=1;
for t=dt:dt:final_t
    for j = 1:3
        [M_trans, RHS_trans] = transientTerm(c_old, dt, alfa+a*k./(1+betta*c));
        M = M_trans-Mdiff+Mbc;
        RHS = RHS_trans+RHSbc;
        c = solvePDE(m,M, RHS);
    end
    c_old = c;
    % visualizeCells(c); drawnow;
    % find the 80% index
    ind_80 = find(c.value(2:end-1)<0.8,1, 'first');
    i=i+1;
    if isempty(ind_80)
        mod_area(i) = 1.0;
    else
        mod_area(i) = 2*x(ind_80)/L;
    end
end
%% visualization
visualizeCells(c);
xlabel('Length [m]'); ylabel('c');

%% getting the 80% modified surface as a function of time
figure
plot(t_mod_wet, mod_area)
