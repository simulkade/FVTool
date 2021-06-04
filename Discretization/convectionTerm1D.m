function M = convectionTerm1D(u)
% This function uses the central difference scheme to discretize a 1D
% convection term in the form \grad (u \phi) where u is a face vactor
% It also returns the x and y parts of the matrix of coefficient.
%
% SYNOPSIS:
%  M = convectionTerm1D(u)
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
G = 1:Nx+2;
DXe = u.domain.cellsize.x(3:end);
DXw = u.domain.cellsize.x(1:end-2);
DXp = u.domain.cellsize.x(2:end-1);

% define the vectors to stores the sparse matrix data
iix = zeros(3*(Nx+2),1);
jjx = zeros(3*(Nx+2),1);
sx = zeros(3*(Nx+2),1);

% reassign the east, west for code readability
ue = u.xvalue(2:Nx+1)./(DXp+DXe);
uw = u.xvalue(1:Nx)./(DXp+DXw);

% calculate the coefficients for the internal cells
AE = reshape(ue,Nx,1);
AW = reshape(-uw,Nx,1);
APx = reshape((ue.*DXe-uw.*DXw)./DXp,Nx,1);

% build the sparse matrix based on the numbering system
rowx_index = reshape(G(2:Nx+1),Nx,1); % main diagonal x
iix(1:3*Nx) = repmat(rowx_index,3,1);
jjx(1:3*Nx) = [reshape(G(1:Nx),Nx,1); ...
		reshape(G(2:Nx+1),Nx,1); reshape(G(3:Nx+2),Nx,1)];
sx(1:3*Nx) = [AW; APx; AE];

% build the sparse matrix
kx = 3*Nx;
M = sparse(iix(1:kx), jjx(1:kx), sx(1:kx), Nx+2, Nx+2);
