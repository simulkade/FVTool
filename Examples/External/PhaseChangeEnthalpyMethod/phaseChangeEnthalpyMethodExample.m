function phaseChangeEnthalpyMethodExample( varargin )
%phaseChangeEnthalpyMethodExample() demonstrates the 2D phase change using 
% the enthalpy method
%   In this example, the heat equation is solved on a 2D initially solid
%   rectangular domain. This file contains three test cases in total. 
%   By setting COMPUTE_FLAG to either 1, 2 or 3, different kinds of phase 
%   change problems for water-ice are solved. Depending on COMPUTE_FLAG, 
%   the visualization is either a plot comparing the numerical solution 
%   with the analytical solution, a full 2D visualization of the 
%   temperature distribution or a full 2D visualization of the phase 
%   distribution within the domain.
%
%   See also StefanProblemAnalyticalSolution.
%
%   Date:   27.09.2015
%   Author: Kai Schueller

close all; % closes all figures
%% Set the flag for the computation:
% COMPUTE_FLAG=1; % Solves the two-phase Stefan problem and plots the
                  % solution together with the analytical solution.
                  % Please note: the analytical solution is for a
                  % semi-infinite domain, while the numerical solution is
                  % for a finite domain. Therefore, the difference between
                  % the solutions increases with increasing time at the
                  % right side of the domain
% COMPUTE_FLAG=2; % Same as COMPUTE_FLAG=1, but 2D visualization without
                  % comparison with analytical solution and Neumann BC
% COMPUTE_FLAG=3; % Same as COMPUTE_FLAG=2, but with convection in the
                  % direction of the heat source (i.e. in negative x) so
                  % that there will be a steady state after some time
COMPUTE_FLAG=3;

%% Add the paths for access to the toolbox functions
addpath('../../../Calculus');
addpath('../../../Visualization');
addpath('../../../Boundary');
addpath('../../../Solvers');
addpath(genpath('../../../Classes'));
addpath('../../../Discretization');
addpath('../../../MeshGeneration');
addpath('../../../Utilities');
addpath('Functions');

%% User Input
% Create the mesh
m=createMesh2D(20,3,2,2);
% visualizeMesh2D(m); % uncomment if mesh should be displayed

% Thermophysical properties of the phase change material (PCM)
rho_L=1000;         % density of the liquid phase [kg/m^3]
rho_S=rho_L;%920;   % density of the solid phase [kg/m^3]
h_melt=333400;      % latent heat of melting [J/kg]
T_m=0;              % melting temperature [degC]
T_L=0.01;           % liquidus temperature [degC]
T_S=-0.01;          % solidus temperature [degC]
k_L=0.6;            % thermal conductivity at liquid phase [W/(mK)]
k_S=2.3;            % thermal conductivity at solid phase [W/(mK)]
c_L=4200;           % Heat capacity liquid Phase [J/(kgK)]
c_S=2000;           % Heat capacity of the solid phase [J/(kgK)]
T_init = -5;

% Create the boundary condition structure
BC = createBC(m);

if COMPUTE_FLAG==1
    % We will later extract a 1D solution, so there is no need for many
    % y-cells
    m=createMesh2D(m.dims(1),2,m.facecenters.x(end),m.facecenters.y(end));
    BC = createBC(m);
    
    u = [0,0];
    T_Solid=-5;     % initial temperature of the solid phase [degC]
    T_Liquid=5;     % initial temperature of the liquid phase [degC]
    analSol.x=0:0.01:2; % x positions for analytical solution
    
    % Assign boundary conditions
    BC.left.a(:) = 0; BC.left.b(:)=1; BC.left.c(:)=T_Liquid;   % Dirichlet for the left boundary 
    BC.right.a(:) = 0; BC.right.b(:)=1; BC.right.c(:)=T_Solid; % Dirichlet for the right boundary
end
if COMPUTE_FLAG==2
    u = [0,0];
    q=-1000; % heat flux [W/m^2]
    % Assign boundary conditions
    BC.left.a(:) = 1; BC.left.b(:)=0; BC.left.c(:)=q;     % Dirichlet for the left boundary 
    BC.right.a(:) = 0; BC.right.b(:)=1; BC.right.c(:)=-5; % Dirichlet for the right boundary 
end
if COMPUTE_FLAG==3
    u = [-0.005,0];      % velocity of the phase change material in x and y
    q=-1000; % heat flux [W/m^2]
    % Assign boundary conditions
    BC.left.a(:) = 1; BC.left.b(:)=0; BC.left.c(:)=q;     % Dirichlet for the left boundary 
    BC.right.a(:) = 0; BC.right.b(:)=1; BC.right.c(:)=-5; % Dirichlet for the right boundary   
end

% Solver settings
dt = 3500; % time step size
final_t = 2000000; % final time
convergence_tolerance=1e-8;
relaxation_coeff=0.7;

%% Set initial values
liquidFraction=createCellVariable(m,0);   % initialize liquid fraction, 
                                          % 0 = initially solid
Temp = createCellVariable(m, T_init, BC); % initial temperature

RHS_liquidFraction=constantSourceTerm(liquidFraction);

%% main loop

% we use constant boundary conditions, so there is no need to change them
% within the main loop
[Mbc, RHSbc] = boundaryCondition2D(BC);

% the convection velocity does also not change
u_face = createFaceVariable(m, u);

for t=dt:dt:final_t
    % The transient term is calculated outside of the liquid fraction
    % update loop for each time step. So the liquid fraction of the last
    % time step is used.
    rho_c_mix=createMixCellVar(m,liquidFraction,rho_L.*c_L,rho_S.*c_S);
    [M_trans, RHS_trans] = transientTerm(Temp, 1, rho_c_mix);
    
    RHS_liquidFraction_old=RHS_liquidFraction;
    
    iterations=0; % set the iteration counter to zero
    
    % set the error higher than tolerance to enter the while-loop for the
    % first time
    error=convergence_tolerance+1;
    
    while error>convergence_tolerance || iterations>10000
        iterations=iterations+1;
        
        % create mixture values
        k_mix=createMixCellVar(m,liquidFraction,k_L,k_S);
        k_mix_face = harmonicMean(k_mix);
        rho_mix=createMixCellVar(m,liquidFraction,rho_L,rho_S);
        rho_mix_face = harmonicMean(rho_mix);
        c_mix=createMixCellVar(m,liquidFraction,c_L,c_S);
        u_rho_mix_face=u_face;
        u_rho_mix_face.xvalue=rho_mix_face.xvalue.*u_face.xvalue;
        u_rho_mix_face.yvalue=rho_mix_face.yvalue.*u_face.yvalue;

        % calculate the matrix M
        Mdiff = dt*diffusionTerm(k_mix_face);
        Mconv =  dt*convectionTerm(u_rho_mix_face);
        M = M_trans-Mdiff+Mbc+Mconv;
        
        % calculate the right hand side vector
        RHS_gammaTerm=constantSourceTerm(rho_mix).*h_melt.*(RHS_liquidFraction_old-RHS_liquidFraction);
        RHS = RHS_trans+RHSbc+RHS_gammaTerm;
        
        Temp_old=Temp; % store the old temperature distribution
        
        % solve the heat equation to obtain the new temperature
        Temp = solvePDE(m,M, RHS);
        
        % calculate maximum error between current and previous temperatures
        error=max(max(abs(Temp.value-Temp_old.value)));
        
        % update the liquid fraction based on the result of heat eq.
        liquidFraction=updateLiquidFraction(liquidFraction,...
                                                Temp,...
                                                c_mix,...
                                                h_melt,...
                                                T_L,...
                                                T_S,...
                                                relaxation_coeff);
        
        RHS_liquidFraction=constantSourceTerm(liquidFraction);
    end
    
    if COMPUTE_FLAG==1
        
        % extract 1D data from solution for the plot
        numSol.T=Temp.value(:,3);
        
        % now interpolate over the ghost cell and first internal cell
        values_temp(1)=(numSol.T(1)+numSol.T(2))/2;
        values_temp(2)=(numSol.T(end)+numSol.T(end-1))/2;
        numSol.T=[values_temp(1); numSol.T(2:end-1); values_temp(2)];
        
        % extracs the x-positions
        numSol.x = [m.facecenters.x(1); m.cellcenters.x; m.facecenters.x(end)];
        
        % calculate analytical solution
        [analSol.T,analSol.InterfacePos] = StefanProblemAnalyticalSolution(...
                                             T_Liquid,...
                                             T_Solid,...
                                             T_m,...
                                             h_melt,...
                                             rho_L,...
                                             k_L,...
                                             k_S,...
                                             c_L,...
                                             c_S,...
                                             analSol.x,...
                                             t);
        
        % plot the analytical and numerical solution
        if t==dt
            h_fig=figure('NumberTitle','Off');
            set(h_fig,'Name',sprintf('Stefan Problem: %d s (%.2f h)',t,t/3600));
            h1=plot(numSol.x,numSol.T,'x','color',[0 0 1]);
            hold on;
            h2=plot(analSol.x,analSol.T,'-','color',[1 0 0]);
            h3=plot([analSol.InterfacePos analSol.InterfacePos],[T_Solid T_Liquid],'--','color',[0 0 0]);
            legend([h1 h2 h3],'Numerical Solution: Temperature','Analytical Solution: Temperature','Analytical Solution: Phase Interface');
            hold off;
        else
            set(h1,'XData',numSol.x,'YData',numSol.T);
            set(h2,'XData',analSol.x,'YData',analSol.T);
            set(h3,'XData',[analSol.InterfacePos analSol.InterfacePos]);
            set(h_fig,'Name',sprintf('Stefan Problem: %d s (%.2f h)',t,t/3600));
        end
        drawnow;
    elseif COMPUTE_FLAG==2 || COMPUTE_FLAG==3
        if COMPUTE_FLAG==2
            visualizeCells(Temp);
            shading interp;
        else
            visualizeCells(liquidFraction);
        end
        
        colormap(jet);
        drawnow; 
    end
end

end

%% Additional Functions

% createMixCellVar() calculates the mixture value of a given property by
% using the liquid fraction of the mixture
function mixCellVar=createMixCellVar(m,liquidFraction,var_L,var_S)
    mixCellVar=liquidFraction.value*var_L+(1-liquidFraction.value)*var_S;
    mixCellVar=createCellVariable(m,mixCellVar(2:end-1,2:end-1));
end

% updateGamma() calculates a new liquid fraction that better fits the
% mixture enthalpy
function helpVar = updateLiquidFraction( liquidFraction_old,T,c,h_melt,T_L,T_S,relaxation_coefficient )
    
    T_corr=liquidFraction_old.value*(T_L-T_S)+T_S;
    
    liquidFraction_new=min(max(liquidFraction_old.value+...
        relaxation_coefficient*c.value./h_melt.*(T.value-T_corr),0),1);
    
    % set liquid fraction of corner cells to zero
    liquidFraction_new([1 end],[1 end])=0;
    helpVar.value=liquidFraction_new;
    helpVar.domain=liquidFraction_old.domain;

end