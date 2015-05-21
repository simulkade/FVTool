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
[Mbc, RHSbc]= boundaryCondition(m,BC);
% define the initial condition
phi0= 0.0; % initial value
phi.Old= createCellVariable(m, phi0, BC); % create initial cells
alfa=createCellVariable(m, 1.0);
phi.value= phi.Old;
% define the time steps
dt= 0.001*L*L/D0; % a proper time step for diffusion process
for i=1:10
    err=100;
    [Mt, RHSt] = transientTerm(m,alfa,dt, phi);
    while err>1e-10
        % calculate the diffusion coefficient
        Dface= harmonicMean(m,D(phi.value));
        % calculate the face value of phi_0
        phi_face= arithmeticMean(m,phi.value);
        % calculate the velocity for convection term
        u= funceval(dD, phi_face).*gradientTerm(m,phi.value);
        % diffusion term
        Mdif= diffusionTerm(m,Dface);
        % convection term
        Mconv= convectionTerm(m,u);
        % divergence term on the RHS
        RHSlin= divergenceTerm(m,u.*phi_face);
        % matrix of coefficients
        M= Mt-Mdif-Mconv+Mbc;
        % RHS vector
        RHS= RHSbc+RHSt-RHSlin;
        % call the linear solver
        phi_new= solvePDE(m, M, RHS);
        % calculate the error
        err= max(abs(phi_new(:)-phi.value(:)));
        % assign the new phi to the old phi
        phi.value=phi_new;
    end
    phi.Old= phi.value;
    visualizeCells(m,phi.value); drawnow;
end