function c_val=tracer_tvd(m, V_dp)
% physical values
mu_val=1e-3; % (Pa.s)
poros=0.2;
perm_val=1.0e-12; % (m^2)
clx=0.05;
cly=0.9;
clz=0.9;
% V_dp=0.6;
if m.dimension==1
    perm=perm_val*ones(m.dims(1),1);
elseif m.dimension==2
    perm= field2d(m.dims(1),m.dims(2),perm_val,V_dp,clx,cly);
elseif
    perm= field3d(m.dims(1),m.dims(2),m.dims(3),perm_val,V_dp,clx,cly,clz);
end
% physical system
u_inj=1.0/(3600*24); % (m/s)
c_init=0.0;
c_inj=1.0;
p_out=100e5; % (Pa)
% create assign values to the domain
k=createCellVariable(m, perm);
labda_face=harmonicMean(k/mu_val);
phi=createCellVariable(m,poros);
% Define the boundaries
BCp = createBC(m); % Neumann BC for pressure
BCc = createBC(m); % Neumann BC for concentration
% change the right boandary to constant pressure (Dirichlet)
BCp.right.a(:)=0.0;
BCp.right.b(:)=1.0;
BCp.right.c(:)=p_out;
% left boundary
BCp.left.a(:)=-labda_face.xvalue(1,:,:);
BCp.left.c(:)=u_inj;
% change the left boundary to constant concentration (Dirichlet)
BCc.left.a(:)=0.0;
BCc.left.b(:)=1.0;
BCc.left.c(:)=1.0;
Mdiffp=diffusionTerm(labda_face);
[Mbcp, RHSp] = boundaryCondition(BCp);
Mp= Mdiffp+Mbcp;
p_val=solvePDE(m, Mp, RHSp);
% estimate a practical time step
n_loop=50;
dt=m.facecenters.x(end)/u_inj/(5*100); % (s)
% find the velocity vector
u=-labda_face.*gradientTerm(p_val); % (m/s)
% find the matrices of coefficients
[Mbcc, RHSbcc] = boundaryCondition(BCc);
% initialize
c_old = createCellVariable(m, c_init, BCc);
c_val = c_old;
FL=fluxLimiter('SUPERBEE');
% start the loop
for i=1:n_loop
    [Mtransc, RHStransc] = transientTerm(c_old, dt, phi);
    for j=1:5 % too lazy to go to full convergence!
        [Mconvc, RHSconv]=convectionTvdTerm(u, c_val, FL);
        Mc=Mconvc+Mbcc+Mtransc;
        RHSc=RHSbcc+RHStransc+RHSconv;
        c_val=solvePDE(m, Mc, RHSc);
    end
    c_old=c_val;
end
end
