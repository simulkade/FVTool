function [M, Mx, My] = convectionTerm2D(u)
% This function uses the central difference scheme to discretize a 2D
% convection term in the form \grad (u \phi) where u is a face vactor
% It also returns the x and y parts of the matrix of coefficient.
%
% SYNOPSIS:
%   [M, Mx, My] = convectionTerm2D(u)
%
% PARAMETERS:
%   u   - FaceVariable  
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
G=reshape(1:(Nx+2)*(Ny+2), Nx+2, Ny+2);
DXe = repmat(u.domain.cellsize.x(3:end), 1, Ny);
DXw = repmat(u.domain.cellsize.x(1:end-2), 1, Ny);
DXp = repmat(u.domain.cellsize.x(2:end-1), 1, Ny);
DYn = repmat(u.domain.cellsize.y(3:end)', Nx, 1);
DYs = repmat(u.domain.cellsize.y(1:end-2)', Nx, 1);
DYp = repmat(u.domain.cellsize.y(2:end-1)', Nx, 1);

% define the vectors to stores the sparse matrix data
iix = zeros(3*(Nx+2)*(Ny+2),1);	iiy = zeros(3*(Nx+2)*(Ny+2),1);
jjx = zeros(3*(Nx+2)*(Ny+2),1);	jjy = zeros(3*(Nx+2)*(Ny+2),1);
sx = zeros(3*(Nx+2)*(Ny+2),1);	sy = zeros(3*(Nx+2)*(Ny+2),1);
mnx = Nx*Ny;	mny = Nx*Ny;

% extract the velocity data
% note: size(ux) = [1:m+1, 1:n] and size(uy) = [1:m, 1:n+1]
ux = u.xvalue;
uy = u.yvalue;

% reassign the east, west, north, and south velocity vectors for the
% code readability
ue = ux(2:Nx+1,:)./(DXp+DXe);		uw = ux(1:Nx,:)./(DXp+DXw);
vn = uy(:,2:Ny+1)./(DYp+DYn);       vs = uy(:,1:Ny)./(DYp+DYs);

% calculate the coefficients for the internal cells
AE = reshape(ue,mnx,1);
AW = reshape(-uw,mnx,1);
AN = reshape(vn,mny,1);
AS = reshape(-vs,mny,1);
APx = reshape((ue.*DXe-uw.*DXw)./DXp,mnx,1);
APy = reshape((vn.*DYn-vs.*DYs)./DYp,mny,1);

% build the sparse matrix based on the numbering system
rowx_index = reshape(G(2:Nx+1,2:Ny+1),mnx,1); % main diagonal x
iix(1:3*mnx) = repmat(rowx_index,3,1);
rowy_index = reshape(G(2:Nx+1,2:Ny+1),mny,1); % main diagonal y
iiy(1:3*mny) = repmat(rowy_index,3,1);
jjx(1:3*mnx) = [reshape(G(1:Nx,2:Ny+1),mnx,1); reshape(G(2:Nx+1,2:Ny+1),mnx,1); reshape(G(3:Nx+2,2:Ny+1),mnx,1)];
jjy(1:3*mny) = [reshape(G(2:Nx+1,1:Ny),mny,1); reshape(G(2:Nx+1,2:Ny+1),mny,1); reshape(G(2:Nx+1,3:Ny+2),mny,1)];
sx(1:3*mnx) = [AW; APx; AE];
sy(1:3*mny) = [AS; APy; AN];

% build the sparse matrix
kx = 3*mnx;
ky = 3*mny;
Mx = sparse(iix(1:kx), jjx(1:kx), sx(1:kx), (Nx+2)*(Ny+2), (Nx+2)*(Ny+2));
My = sparse(iiy(1:ky), jjy(1:ky), sy(1:ky), (Nx+2)*(Ny+2), (Nx+2)*(Ny+2));
M = Mx + My;
