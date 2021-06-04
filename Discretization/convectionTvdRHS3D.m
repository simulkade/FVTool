function [RHS, RHSx, RHSy, RHSz] = ...
    convectionTvdRHS3D(u, phi, FL)
% This function uses the TVD scheme to discretize a 3D
% convection term in the form \grad (u \phi) where u is a face vactor
% It also returns the x and y parts of the matrix of coefficient.
%
% SYNOPSIS:
%   [M, RHS, Mx, My, Mz, RHSx, RHSy, RHSz] = ...
%    convectionTvdTerm3D(u, phi, FL)
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
Nx = u.domain.dims(1);
Ny = u.domain.dims(2);
Nz = u.domain.dims(3);
G=reshape((1:(Nx+2)*(Ny+2)*(Nz+2)), Nx+2, Ny+2, Nz+2);
DXp = repmat(u.domain.cellsize.x(2:end-1), 1, Ny, Nz);
DYp = repmat(u.domain.cellsize.y(2:end-1)', Nx, 1, Nz);
DZ = zeros(1,1,Nz+2);
DZ(1,1,:) = u.domain.cellsize.z;
DZp=repmat(DZ(1,1,2:end-1), Nx, Ny, 1);
dx=repmat(0.5*(u.domain.cellsize.x(1:end-1)+u.domain.cellsize.x(2:end)), 1, Ny, Nz);
dy=repmat(0.5*(u.domain.cellsize.y(1:end-1)+u.domain.cellsize.y(2:end))', Nx, 1, Nz);
dz=zeros(1, 1, Nz+1);
dz(1,1,:)=0.5*(u.domain.cellsize.z(1:end-1)+u.domain.cellsize.z(2:end));
dz=repmat(dz, Nx, Ny, 1);
psiX_p = zeros(Nx+1,Ny,Nz);
psiX_m = zeros(Nx+1,Ny,Nz);
psiY_p = zeros(Nx,Ny+1,Nz);
psiY_m = zeros(Nx,Ny+1,Nz);
psiZ_p = zeros(Nx,Ny,Nz+1);
psiZ_m = zeros(Nx,Ny,Nz+1);

% define the vectors to stores the sparse matrix data
mnx = Nx*Ny*Nz;	mny = Nx*Ny*Nz;   mnz = Nx*Ny*Nz;

% extract the velocity data
% note: size(ux) = [1:m+1, 1:n] and size(uy) = [1:m, 1:n+1]
ux = u.xvalue;
uy = u.yvalue;
uz = u.zvalue;

% calculate the upstream to downstream gradient ratios for u>0 (+ ratio)
% x direction
dphiX_p = (phi.value(2:Nx+2, 2:Ny+1, 2:Nz+1)-phi.value(1:Nx+1, 2:Ny+1, 2:Nz+1))./dx;
rX_p = dphiX_p(1:end-1,:,:)./fsign(dphiX_p(2:end,:,:));
psiX_p(2:Nx+1,:,:) = 0.5*FL(rX_p).* ...
    (phi.value(3:Nx+2,2:Ny+1,2:Nz+1)-phi.value(2:Nx+1,2:Ny+1,2:Nz+1));
psiX_p(1,:,:) = 0; % left boundary
% y direction
dphiY_p = (phi.value(2:Nx+1, 2:Ny+2, 2:Nz+1)-phi.value(2:Nx+1, 1:Ny+1, 2:Nz+1))./dy;
rY_p = dphiY_p(:,1:end-1,:)./fsign(dphiY_p(:,2:end,:));
psiY_p(:,2:Ny+1,:) = 0.5*FL(rY_p).* ...
    (phi.value(2:Nx+1,3:Ny+2,2:Nz+1)-phi.value(2:Nx+1, 2:Ny+1,2:Nz+1));
psiY_p(:,1,:) = 0; % Bottom boundary
% z direction
dphiZ_p = (phi.value(2:Nx+1, 2:Ny+1, 2:Nz+2)-phi.value(2:Nx+1, 2:Ny+1, 1:Nz+1))./dz;
rZ_p = dphiZ_p(:,:,1:end-1)./fsign(dphiZ_p(:,:,2:end));
psiZ_p(:,:,2:Nz+1) = 0.5*FL(rZ_p).* ...
    (phi.value(2:Nx+1,2:Ny+1,3:Nz+2)-phi.value(2:Nx+1,2:Ny+1,2:Nz+1));
psiZ_p(:,:,1) = 0; % Back boundary

% calculate the upstream to downstream gradient ratios for u<0 (- ratio)
% x direction
rX_m = dphiX_p(2:end,:,:)./fsign(dphiX_p(1:end-1,:,:));
psiX_m(1:Nx,:,:) = 0.5*FL(rX_m).* ...
    (phi.value(1:Nx, 2:Ny+1, 2:Nz+1)-phi.value(2:Nx+1, 2:Ny+1, 2:Nz+1));
psiX_m(Nx+1,:,:) = 0; % right boundary
% y direction
rY_m = dphiY_p(:,2:end,:)./fsign(dphiY_p(:,1:end-1,:));
psiY_m(:,1:Ny,:) = 0.5*FL(rY_m).* ...
    (phi.value(2:Nx+1,1:Ny,2:Nz+1)-phi.value(2:Nx+1,2:Ny+1,2:Nz+1));
psiY_m(:,Ny+1,:) = 0; % top boundary
% z direction
rZ_m = dphiZ_p(:,:,2:end)./fsign(dphiZ_p(:,:,1:end-1));
psiZ_m(:,:,1:Nz) = 0.5*FL(rZ_m).* ...
    (phi.value(2:Nx+1,2:Ny+1,1:Nz)-phi.value(2:Nx+1,2:Ny+1,2:Nz+1));
psiZ_m(:,:,Nz+1) = 0; % front boundary
% reassign the east, west, north, and south velocity vectors for the
% code readability
ue = ux(2:Nx+1,:,:);		uw = ux(1:Nx,:,:);
vn = uy(:,2:Ny+1,:);     vs = uy(:,1:Ny,:);
wf = uz(:,:,2:Nz+1);     wb = uz(:,:,1:Nz);

% find the velocity direction for the upwind scheme
ue_min = min(ue,0);	ue_max = max(ue,0);
uw_min = min(uw,0);	uw_max = max(uw,0);
vn_min = min(vn,0);	vn_max = max(vn,0);
vs_min = min(vs,0);	vs_max = max(vs,0);
wf_min = min(wf,0);	wf_max = max(wf,0);
wb_min = min(wb,0);	wb_max = max(wb,0);

% build the sparse matrix based on the numbering system
rowx_index = reshape(G(2:Nx+1,2:Ny+1,2:Nz+1),mnx,1); % main diagonal x
rowy_index = reshape(G(2:Nx+1,2:Ny+1,2:Nz+1),mny,1); % main diagonal y
rowz_index = reshape(G(2:Nx+1,2:Ny+1,2:Nz+1),mnz,1); % main diagonal z

% calculate the TVD correction term
div_x = -(1./DXp).*((ue_max.*psiX_p(2:Nx+1,:,:)+ue_min.*psiX_m(2:Nx+1,:,:))- ...
              (uw_max.*psiX_p(1:Nx,:,:)+uw_min.*psiX_m(1:Nx,:,:)));
div_y = -(1./DYp).*((vn_max.*psiY_p(:,2:Ny+1,:)+vn_min.*psiY_m(:,2:Ny+1,:))- ...
              (vs_max.*psiY_p(:,1:Ny,:)+vs_min.*psiY_m(:,1:Ny,:)));
div_z = -(1./DZp).*((wf_max.*psiZ_p(:,:,2:Nz+1)+wf_min.*psiZ_m(:,:,2:Nz+1))- ...
              (wb_max.*psiZ_p(:,:,1:Nz)+wb_min.*psiZ_m(:,:,1:Nz)));

% define the RHS Vector
RHS = zeros((Nx+2)*(Ny+2)*(Nz+2),1);
RHSx = zeros((Nx+2)*(Ny+2)*(Nz+2),1);
RHSy = zeros((Nx+2)*(Ny+2)*(Nz+2),1);
RHSz = zeros((Nx+2)*(Ny+2)*(Nz+2),1);

% assign the values of the RHS vector
row_index = rowx_index;
RHS(row_index) = reshape(div_x+div_y+div_z,Nx*Ny*Nz,1);
RHSx(rowx_index) = reshape(div_x,Nx*Ny*Nz,1);
RHSy(rowy_index) = reshape(div_y,Nx*Ny*Nz,1);
RHSz(rowz_index) = reshape(div_z,Nx*Ny*Nz,1);

end

function phi_out = fsign(phi_in)
% This function checks the value of phi_in and assigns an eps value to the
% elements that are less than or equal to zero, while keeping the signs of
% the nonzero elements
    phi_out = (abs(phi_in)>=eps).*phi_in+eps*(phi_in==0)+eps*(abs(phi_in)<eps).*sign(phi_in);
end
