% Nonlinear diffusion equation: a tutorial
% see http://fvt.simulkade.com/posts/2015-04-06-solving-nonlinear-pdes-with-fvm.html
L= 1.0; % domain length
Nx= 100; % number of cells
m= createMesh1D(Nx, L); % create a 1D mesh
D0= 1.0; % diffusion coefficient constant
% Define the diffusion coefficientand its derivative
D=@(phi)(D0*(1.0+phi.^2));
dD=@(phi)(2.0*D0*phi);
% create boundary condition
BC = createBC(m);
BC.left.a(:)=0.0;
BC.left.b(:)=1.0;
BC.left.c(:)=5.0;
BC.right.a(:)=0.0;
BC.right.b(:)=1.0;
BC.right.c(:)=0.0;
[Mbc, RHSbc]= boundaryCondition(BC);
% define the initial condition
phi0= 0.0; % initial value
phi_old= createCellVariable(m, phi0, BC); % create initial cells
alfa=createCellVariable(m, 1.0);
phi= phi_old;
% define the time steps
dt= 0.001*L*L/D0; % a proper time step for diffusion process
for i=1:10
    err=100;
    [Mt, RHSt] = transientTerm(phi_old, dt, alfa);
    while err>1e-10
        % calculate the diffusion coefficient
        Dface= harmonicMean(D(phi));
        % calculate the face value of phi_0
        phi_face= arithmeticMean(phi);
        % calculate the velocity for convection term
        u= funceval(dD, phi_face).*gradientTerm(phi);
        % diffusion term
        Mdif= diffusionTerm(Dface);
        % convection term
        Mconv= convectionTerm(u);
        % divergence term on the RHS
        RHSlin= divergenceTerm(u.*phi_face);
        % matrix of coefficients
        M= Mt-Mdif-Mconv+Mbc;
        % RHS vector
        RHS= RHSbc+RHSt-RHSlin;
        % call the linear solver
        phi_new= solvePDE(m, M, RHS);
        % calculate the error
        err= max(abs(phi_new.value(:)-phi.value(:)));
        % assign the new phi to the old phi
        phi=phi_new;
    end
    phi_old= phi;
    visualizeCells(phi); drawnow;
end
