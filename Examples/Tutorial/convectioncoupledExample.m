% 1D convection coupled example with source terms
% diffusion coefficients are defined but not used
% three equations are coupled
% scheme can be
clc; clear;
% define a 1D domain and mesh
W = 1;
Nx = 300;
mesh1 = createMesh1D(Nx, W);
x = mesh1.cellcenters.x;
% define the boundaries
BC = createBC(mesh1); % all Neumann
BC.left.a(:)=0; BC.left.b(:)=1; BC.left.c(:)=1;
BC.right.a(:) =0; BC.right.b(:) =1; BC.right.c(:) =0;
% Initial values
phi1_old = createCellVariable(mesh1, 0.0, BC);
phi2_old = createCellVariable(mesh1, 0.0, BC);
phi3_old = createCellVariable(mesh1, 0.0, BC);
% initial guess for phi
phi1 = phi1_old;
phi2 = phi2_old;
phi3 = phi3_old;
% velocity field
V= [6.9444; 3.4722; 2.3148];
uf1 = createFaceVariable(mesh1, V(1));
uf2 = createFaceVariable(mesh1, V(2));
uf3 = createFaceVariable(mesh1, V(3));
% diffusion field
D = 1e-2;
Df = createFaceVariable(mesh1, D);
% transient term coefficient
alfa = createCellVariable(mesh1,1.0);
% upwind convection term
Mconv1 = convectionUpwindTerm(uf1);
Mconv2 = convectionUpwindTerm(uf2);
Mconv3 = convectionUpwindTerm(uf3);
% define the BC term
[Mbc, RHSbc] = boundaryCondition(BC);
% choose a flux limiter
FL = fluxLimiter('Superbee');
% source terms
M=[-0.1200, 0.1200, 0; 0.0579, -0.1779, 0.1200; 0, 0.0772, -0.0772];
Ms=linearSourceTerm(alfa);
% solver
dt = 0.0005; % time step
final_t = W/max(V);
t = 0;
while t<final_t
    t = t+dt;
    % inner loop for TVD scheme
    [Mt1, RHSt1] = transientTerm(phi1_old, dt, alfa);
    [Mt2, RHSt2] = transientTerm(phi2_old, dt, alfa);
    [Mt3, RHSt3] = transientTerm(phi3_old, dt, alfa);
    Muw=[Mt1+Mconv1-M(1,1)*Ms+Mbc -M(1,2)*Ms -M(1,3)*Ms
         -M(2,1)*Ms Mt2+Mconv2-M(2,2)*Ms+Mbc  -M(2,3)*Ms
         -M(3,1)*Ms -M(3,2)*Ms Mt2+Mconv3-M(3,3)*Ms+Mbc];
    RHS=[RHSt1+RHSbc; RHSt2+RHSbc; RHSt3+RHSbc];
    X=Muw\RHS;
    phi1.value=X(1:Nx+2);
    phi2.value=X((Nx+2+1):(2*(Nx+2)));
    phi3.value=X((2*(Nx+2)+1):end);
    phi1_old = phi1;
    phi2_old = phi2;
    phi3_old = phi3;
    figure(1);plot(x, phi1.value(2:Nx+1), x, phi2.value(2:Nx+1), '-o', x, ...
        phi3.value(2:Nx+1)); drawnow;

end
