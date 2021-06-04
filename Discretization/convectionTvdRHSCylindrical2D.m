function [RHS, RHSx, RHSy] = ...
    convectionTvdRHSCylindrical2D(u, phi, FL)
% This function uses the upwind scheme to discretize a 1D
% convection term in the form \grad (u \phi) where u is a face vactor
% It also returns the x and y parts of the matrix of coefficient.
%
% SYNOPSIS:
%   [M, RHS, Mx, My, RHSx, RHSy] = ...
%    convectionTvdTermCylindrical2D(u, phi, FL)
%
% PARAMETERS:
%
%
% RETURNS:
%
%
% EXAMPLE:
%
% SEE ALSO:
%

% extract data from the mesh structure
Nr = u.domain.dims(1);
Nz = u.domain.dims(2);
G=reshape(1:(Nr+2)*(Nz+2), Nr+2, Nz+2);
DRp = repmat(u.domain.cellsize.x(2:end-1), 1, Nz);
DZp = repmat(u.domain.cellsize.y(2:end-1)', Nr, 1);
rp = repmat(u.domain.cellcenters.x, 1, Nz);
rf = repmat(u.domain.facecenters.x, 1, Nz);
dr=repmat(0.5*(u.domain.cellsize.x(1:end-1)+u.domain.cellsize.x(2:end)), 1, Nz);
dz=repmat(0.5*(u.domain.cellsize.y(1:end-1)+u.domain.cellsize.y(2:end))', Nr, 1);
psiX_p = zeros(Nr+1,Nz);
psiX_m = zeros(Nr+1,Nz);
psiY_p = zeros(Nr,Nz+1);
psiY_m = zeros(Nr,Nz+1);

% define the vectors to stores the sparse matrix data
mNr = Nr*Nz;	mNz = Nr*Nz;

% extract the velocity data
% note: size(ux) = [1:m+1, 1:n] and size(uy) = [1:m, 1:n+1]
ux = u.xvalue;
uy = u.yvalue;

% calculate the upstream to downstream gradient ratios for u>0 (+ ratio)
% x direction
dphiX_p = (phi.value(2:Nr+2, 2:Nz+1)-phi.value(1:Nr+1, 2:Nz+1))./dr;
rX_p = dphiX_p(1:end-1,:)./fsign(dphiX_p(2:end,:));
psiX_p(2:Nr+1,:) = 0.5*FL(rX_p).*(phi.value(3:Nr+2,2:Nz+1)-phi.value(2:Nr+1, 2:Nz+1));
psiX_p(1, :) = 0; % left boundary will be handled in the main matrix
% y direction
dphiY_p = (phi.value(2:Nr+1, 2:Nz+2)-phi.value(2:Nr+1, 1:Nz+1))./dz;
rY_p = dphiY_p(:,1:end-1)./fsign(dphiY_p(:,2:end));
psiY_p(:,2:Nz+1) = 0.5*FL(rY_p).*(phi.value(2:Nr+1,3:Nz+2)-phi.value(2:Nr+1, 2:Nz+1));
psiY_p(:,1) = 0; % Bottom boundary will be handled in the main matrix

% calculate the upstream to downstream gradient ratios for u<0 (- ratio)
% x direction
rX_m = dphiX_p(2:end,:)./fsign(dphiX_p(1:end-1,:));
psiX_m(1:Nr,:) = 0.5*FL(rX_m).*(phi.value(1:Nr, 2:Nz+1)-phi.value(2:Nr+1, 2:Nz+1));
psiX_m(Nr+1,:) = 0; % right boundary
% y direction
rY_m = dphiY_p(:,2:end)./fsign(dphiY_p(:,1:end-1));
psiY_m(:,1:Nz) = 0.5*FL(rY_m).*(phi.value(2:Nr+1, 1:Nz)-phi.value(2:Nr+1, 2:Nz+1));
psiY_m(:, Nz+1) = 0; % top boundary will be handled in the main matrix

% reassign the east, west, north, and south velocity vectors for the
% code readability
ue = ux(2:Nr+1,:);		uw = ux(1:Nr,:);
vn = uy(:,2:Nz+1);       vs = uy(:,1:Nz);
re = rf(2:Nr+1,:);         rw = rf(1:Nr,:);

% find the velocity direction for the upwind scheme
ue_min = min(ue,0);	ue_max = max(ue,0);
uw_min = min(uw,0);	uw_max = max(uw,0);
vn_min = min(vn,0);	vn_max = max(vn,0);
vs_min = min(vs,0);	vs_max = max(vs,0);

% build the sparse matrix based on the numbering system
rowx_index = reshape(G(2:Nr+1,2:Nz+1),mNr,1); % main diagonal x
rowy_index = reshape(G(2:Nr+1,2:Nz+1),mNz,1); % main diagonal y

% calculate the TVD correction term
div_x = -(1./(DRp.*rp)).*(re.*(ue_max.*psiX_p(2:Nr+1,:)+ue_min.*psiX_m(2:Nr+1,:))- ...
              rw.*(uw_max.*psiX_p(1:Nr,:)+uw_min.*psiX_m(1:Nr,:)));
div_y = -(1./DZp).*((vn_max.*psiY_p(:,2:Nz+1)+vn_min.*psiY_m(:,2:Nz+1))- ...
              (vs_max.*psiY_p(:,1:Nz)+vs_min.*psiY_m(:,1:Nz)));
% define the RHS Vector
RHS = zeros((Nr+2)*(Nz+2),1);
RHSx = zeros((Nr+2)*(Nz+2),1);
RHSy = zeros((Nr+2)*(Nz+2),1);

% assign the values of the RHS vector
RHS(rowx_index) = reshape(div_x+div_y,Nr*Nz,1);
RHSx(rowx_index) = reshape(div_x,Nr*Nz,1);
RHSy(rowy_index) = reshape(div_y,Nr*Nz,1);

end

function phi_out = fsign(phi_in)
% This function checks the value of phi_in and assigns an eps value to the
% elements that are less than or equal to zero, while keeping the signs of
% the nonzero elements
    phi_out = (abs(phi_in)>=eps).*phi_in+eps*(phi_in==0)+eps*(abs(phi_in)<eps).*sign(phi_in);
end
