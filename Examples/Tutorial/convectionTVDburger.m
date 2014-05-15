% BURGER's equation; there are still a few strange outputs when I use TVD
% but in general works just fine, wait for a blog post
clc; clear;
% define a 1D mesh
W = 1;
Nx = 600;
mesh1 = buildMesh1D(Nx, W);
x = mesh1.cellcenters.x;
% define the boundaries
BC = createBC(mesh1); % all Neumann
BC.left.periodic=1;
BC.right.periodic =1;
% Initial values
phi.Old = createCellVariable(mesh1, 0, BC);
phi.Old(2:floor(end/2)) = sin(x(1:floor(end/2))*2*pi());
% phi.Old(180:400)= sin(x(180:400)*10*pi());
phi.value = phi.Old;
phiuw = phi;
phiinit=phi.Old;
% velocity field
u = 1;
uf = createFaceVariable(mesh1, u);
% uf.xvalue = linspace(0.4,0.2, Nx+1)';
% transient term coefficient
alfa = ones(Nx,1);
% matrix of coefficients
Mconvuw = convectionUpwindTerm1D(mesh1, uf);
% define the BC term
[Mbc, RHSbc] = boundaryCondition(mesh1, BC);
% choose a flux limiter
FL = fluxLimiter('MUSCL');
phi_face = arithmeticMean(mesh1, phi.value);
% solver
dt = 0.001;
t = 0;
for i = 1:100000
    t = t+dt;
    % define the transient term
    for j = 1:10
        [Mt, RHSt] = transientTerm(mesh1, alfa, dt, phi);
        phi_face = tvdMean(mesh1, phi.value, uf, FL);
        uf = phi_face;
        [Mconv, RHSconv] = convectionTvdTerm(mesh1, uf, phi.value, FL);
        RHSdiv = divergenceTerm(mesh1, 0.5*phi_face.*phi_face);
    %     Mconv = convectionUpwindTerm1D(mesh1, uf);
        M = Mconv+Mt+Mbc;
        RHS = RHSt+RHSbc+RHSconv+RHSdiv;
        PHI = M\RHS;
        phi.value = reshape(PHI, Nx+2, 1);
    end
    for j = 1:5
        [Mtuw, RHStuw] = transientTerm(mesh1, alfa, dt, phiuw);
        phi_face = upwindMean(mesh1, uf, phiuw.value);
        uf = phi_face;
        Mconvuw = convectionUpwindTerm1D(mesh1, uf);
        RHSdiv = divergenceTerm(mesh1, 0.5*phi_face.*phi_face);
        Muw = Mconvuw+Mtuw+Mbc;
        RHSuw = RHStuw+RHSbc+RHSdiv;
        PHI = Muw\RHSuw;
        phiuw.value = reshape(PHI, Nx+2, 1);
    end
    phiuw.Old = phiuw.value;
    phi.Old = phi.value;
    figure(1);plot(x, phiinit(2:Nx+1), x, phi.value(2:Nx+1), '-o', x, phiuw.value(2:Nx+1)); drawnow;
end
