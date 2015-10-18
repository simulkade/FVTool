function [ U_cellVar,V_cellVar,U_aP_sumMat,V_aP_sumMat ] = MomentumEq(veloRelax,m, pGradMat,mu_faceVar,faceVelocity,U_old,V_old,U_BC,V_BC,rho,ALGORITHM_TYPE )
%MOMENTUMEQ Summary of this function goes here
%   Detailed explanation goes here

% Discretization of the diffusive term
Mdiff = diffusionTerm(mu_faceVar);

% Discretization of the convective term
% could be improved by using the power law scheme
Mconv = rho*convectionUpwindTerm(faceVelocity);

% Boundary Conditions
[U_Mbc, U_RHSbc] = boundaryCondition2D(U_BC);
[V_Mbc, V_RHSbc] = boundaryCondition2D(V_BC);

% build Matrices for U and V
U_M = -Mdiff+U_Mbc+Mconv;
V_M = -Mdiff+V_Mbc+Mconv;

% extract aP
% it looks like U_ap and V_ap are the same, but i'm not sure if this is the
% case for all BC's. So we could also just use U_ap, if speedup is required
U_aP_sumMat=calcApSum(m,U_M,ALGORITHM_TYPE);
V_aP_sumMat=calcApSum(m,V_M,ALGORITHM_TYPE);

% check if they are the same
%if find(U_aP_sum-V_aP_sum)
    %errordlg('Fehler');
%end

% create separate cell variables for gradient in x and y
% directions, because that's what constantSourceTerm expects
pGradxCellVar=createCellVariable(m,0);
pGradyCellVar=pGradxCellVar;
pGradxCellVar.value(2:end-1,2:end-1)=pGradMat.xvalue;
pGradyCellVar.value(2:end-1,2:end-1)=pGradMat.yvalue;

% build RHS vector for U and V
U_RHS = U_RHSbc-constantSourceTerm(pGradxCellVar)+(1-veloRelax)./veloRelax.*diag(U_M).*reshape(U_old.value,[],1);
V_RHS = V_RHSbc-constantSourceTerm(pGradyCellVar)+(1-veloRelax)./veloRelax.*diag(V_M).*reshape(V_old.value,[],1);

U_M(linspace(1,numel(U_M),length(U_M)))  = diag(U_M)/veloRelax;
V_M(linspace(1,numel(V_M),length(V_M)))  = diag(V_M)/veloRelax;

U_cellVar = solvePDE(m,U_M, U_RHS);
V_cellVar = solvePDE(m,V_M, V_RHS);

end

function ApSumMat=calcApSum(m,M,ALGORITHM_TYPE)
    
    Nx = m.dims(1);
    Ny = m.dims(2);
    epsilon=0; % this can be set to a small value to prevent division by zero

    switch lower(ALGORITHM_TYPE)
        case 'simple'
            ApSum=diag(M);
        case 'simplec' % not working yet
            ApSum=(2*diag(M)+full(sum(M,2)))*(m.cellsize.x(1)*m.cellsize.y(1))+epsilon;
    end
    
    % it is necessary to multiply the coefficients with the volume to be in
    % accordance with the literature, because with FVTool the coefficients
    % are already divided by the volume
    volume=repmat(m.cellsize.x,1,Ny+2).*repmat(m.cellsize.y',Nx+2,1);
    ApSumMat.value = reshape(ApSum,Nx+2,Ny+2).*volume+epsilon;
    ApSumMat.domain=m;  
    
end

