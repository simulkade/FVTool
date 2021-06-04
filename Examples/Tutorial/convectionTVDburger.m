% BURGER's equation; there are still a few strange outputs when I use TVD
% but in general works just fine, wait for a blog post
% The solution does not look bad and kind of satisfying to look at.
% increate the number of iteration if you plan to bore yourself to death.
% Written by Ali A. Eftekhari
% Last checked: June 2021
clc;
% define a 1D mesh
W = 1;
Nx = 600;
mesh1 = createMesh1D(Nx, W);
x = mesh1.cellcenters.x;
% define the boundaries
BC = createBC(mesh1); % all Neumann
BC.left.periodic=1;
BC.right.periodic =1;
% Initial values
phi_old = createCellVariable(mesh1, 0, BC);
phi_old.value(2:floor(end/2)) = sin(x(1:floor(end/2))*2*pi());
% phi.Old(180:400)= sin(x(180:400)*10*pi());
phi = phi_old;
phiuw = phi;
phiinit=phi_old;
phiuw_old=phi;
% velocity field
u = 1;
uf = createFaceVariable(mesh1, u);
% uf.xvalue = linspace(0.4,0.2, Nx+1)';
% transient term coefficient
alfa = createCellVariable(mesh1,1);
% matrix of coefficients
Mconvuw = convectionUpwindTerm1D(uf);
% define the BC term
[Mbc, RHSbc] = boundaryCondition(BC);
% choose a flux limiter
FL = fluxLimiter('Superbee');
phi_face = linearMean(phi);
% solver
dt = 0.001;
t = 0;
for i = 1:1000
    t = t+dt;
    % define the transient term
    [Mt, RHSt] = transientTerm(phi_old, dt, alfa);
    for j = 1:10
        phi_face = tvdMean(phi, uf, FL);
        uf = phi_face;
        [Mconv, RHSconv] = convectionTvdTerm(uf, phi, FL);
        RHSdiv = divergenceTerm(0.5*phi_face.*phi_face);
    %     Mconv = convectionUpwindTerm1D(mesh1, uf);
        M = Mconv+Mt+Mbc;
        RHS = RHSt+RHSbc+RHSconv+RHSdiv;
        phi = solvePDE(mesh1, M, RHS);
    end
    [Mtuw, RHStuw] = transientTerm(phiuw_old, dt, alfa);
    for j = 1:5
        phi_face = upwindMean(phiuw, uf);
        uf = phi_face;
        Mconvuw = convectionUpwindTerm1D(uf);
        RHSdiv = divergenceTerm(0.5*phi_face.*phi_face);
        Muw = Mconvuw+Mtuw+Mbc;
        RHSuw = RHStuw+RHSbc+RHSdiv;
        phiuw = solvePDE(mesh1, Muw, RHSuw);
    end
    phiuw_old = phiuw;
    phi_old = phi;
    figure(1);plot(x, phiinit.value(2:Nx+1), x, phi.value(2:Nx+1), '-o', x, phiuw.value(2:Nx+1)); drawnow;
end
