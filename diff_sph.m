clc; clear;
% define parameters
c0 = 1;

% define a 1D mesh
L = 1;
Nr = 1000;
mesh1 = createMeshSpherical1D(Nr, L);
x = mesh1.cellcenters.x;
xf = mesh1.facecenters.x;

% define the boundaries
BC = createBC(mesh1); % all Neumann
BC.right.a=0;
BC.right.b=1;
BC.right.c=0;
[Mbc, RHSbc] = boundaryCondition(BC);

% define diffusion term
D = createCellVariable(mesh1,1);
Dave = harmonicMean(D);
Mdiff = diffusionTermSpherical1D(Dave);

% initial values
c_old = createCellVariable(mesh1, c0);
c = c_old;

% solver
dt = 1E-3;
t = dt;
N = 30;

err = 0;
tvec = dt;

it = 1;

while t < 1
    % define the transient term
    [Mt, RHSt] = transientTerm(c_old, dt);
    
    % assemble the matrix and vector
    M = Mt-Mdiff+Mbc;
    RHS = RHSt+RHSbc;
    
    % solve PDE
    c = solvePDE(mesh1,M,RHS);
    c_old = c;
    
    % calculate analytic solution
    f = zeros(Nr,1);
    for n=1:N
        f = f + 2*(-1)^(n+1)/(n*pi)*exp(-(n*pi)^2*t)*sin(n*pi*x)./x;
    end
    
    tvec = [tvec; t];
    err = [err; trapz(x,(f-c.value(2:Nr+1)).^2)];
    
    if mod(it,10) == 0
    figure(1);
    plot(x, f, 'b', [0; x; 1], c.value,'ro'); 
    title(sprintf('t=%.2e',t));
    xlabel('r')
    ylabel('c')
    legend('analytical','numerical')
    ylim([0,1]); 
    drawnow;
    end
    t = t+dt;
    it = it + 1;
end

figure(2)
plot(tvec,sqrt(err))
xlabel('t')
ylabel('absolute error')