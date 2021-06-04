function [RHSdiv, RHSdivx, RHSdivy, RHSdivz] = divergenceTermCylindrical3D(F)
% This function calculates the divergence of a field using its face
% average flux vector F, which is a face vector
%
% SYNOPSIS:
%   [RHSdiv, RHSdivx, RHSdivy, RHSdivz] = divergenceTermCylindrical3D(F)
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
Nx = F.domain.dims(1);
Ntheta = F.domain.dims(2);
Nz = F.domain.dims(3);
G=reshape((1:(Nx+2)*(Ntheta+2)*(Nz+2)), Nx+2, Ntheta+2, Nz+2);
dx = repmat(F.domain.cellsize.x(2:end-1), 1, Ntheta, Nz);
dy = repmat(F.domain.cellsize.y(2:end-1)', Nx, 1, Nz);
DZ = zeros(1,1,Nz+2);
DZ(1,1,:) = F.domain.cellsize.z;
dz=repmat(DZ(1,1,2:end-1), Nx, Ntheta, 1);
rp = repmat(F.domain.cellcenters.x, 1, Ntheta, Nz);

% define the vector of cell index
row_index = reshape(G(2:Nx+1,2:Ntheta+1, 2:Nz+1),Nx*Ntheta*Nz,1); % main diagonal (only internal cells)

% calculate the flux vector in x and y direction
% note: size(Fx) = [1:m+1, 1:n] and size(Fy) = [1:m, 1:n+1]
Fx = F.xvalue;
Fy = F.yvalue;
Fz = F.zvalue;

% reassign the east, west, north, and south flux vectors for the
% code readability
Fe = Fx(2:Nx+1,:,:);		Fw = Fx(1:Nx,:,:);
Fn = Fy(:,2:Ntheta+1,:);     Fs = Fy(:,1:Ntheta,:);
Ff = Fz(:,:,2:Nz+1);     Fb = Fz(:,:,1:Nz);
% compute the divergence
div_x = (Fe-Fw)./dx;
div_y = (Fn-Fs)./(dy.*rp);
div_z = (Ff-Fb)./dz;

% define the RHS Vector
RHSdiv = zeros((Nx+2)*(Ntheta+2)*(Nz+2),1);
RHSdivx = zeros((Nx+2)*(Ntheta+2)*(Nz+2),1);
RHSdivy = zeros((Nx+2)*(Ntheta+2)*(Nz+2),1);
RHSdivz = zeros((Nx+2)*(Ntheta+2)*(Nz+2),1);

% assign the values of the RHS vector
RHSdiv(row_index) = reshape(div_x+div_y+div_z,Nx*Ntheta*Nz,1);
RHSdivx(row_index) = reshape(div_x,Nx*Ntheta*Nz,1);
RHSdivy(row_index) = reshape(div_y,Nx*Ntheta*Nz,1);
RHSdivz(row_index) = reshape(div_z,Nx*Ntheta*Nz,1);
