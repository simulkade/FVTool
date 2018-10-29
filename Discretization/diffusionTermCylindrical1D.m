function M = diffusionTermCylindrical1D(D)
% This function uses the central difference scheme to discretize a 1D
% diffusion term in the form \grad . (D \grad \phi) where u is a face vactor
% It also returns the x and y parts of the matrix of coefficient.
%
% SYNOPSIS:
%
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

% Copyright (c) 2012-2016 Ali Akbar Eftekhari
% See the license file

% extract data from the mesh structure
Nx = D.domain.dims(1);
G = 1:Nx+2;
DX = D.domain.cellsize.x;
dx = 0.5*(DX(1:end-1)+DX(2:end));
rp = D.domain.cellcenters.x;
rf = D.domain.facecenters.x;

% define the vectors to store the sparse matrix data
iix = zeros(3*(Nx+2),1);
jjx = zeros(3*(Nx+2),1);
sx = zeros(3*(Nx+2),1);

% extract the velocity data
% note: size(Dx) = [1:m+1, 1:n] and size(Dy) = [1:m, 1:n+1]
Dx = D.xvalue;

% reassign the east, west, north, and south velocity vectors for the
% code readability
De = rf(2:Nx+1).*Dx(2:Nx+1)./(rp.*dx(2:Nx+1).*DX(2:Nx+1));
Dw = rf(1:Nx).*Dx(1:Nx)./(rp.*dx(1:Nx).*DX(2:Nx+1));

% calculate the coefficients for the internal cells
AE = reshape(De,Nx,1);
AW = reshape(Dw,Nx,1);
APx = reshape(-(De+Dw),Nx,1);

% build the sparse matrix based on the numbering system
rowx_index = reshape(G(2:Nx+1),Nx,1); % main diagonal x
iix(1:3*Nx) = repmat(rowx_index,3,1);
jjx(1:3*Nx) = [reshape(G(1:Nx),Nx,1); reshape(G(2:Nx+1),Nx,1); reshape(G(3:Nx+2),Nx,1)];
sx(1:3*Nx) = [AW; APx; AE];

% build the sparse matrix
kx = 3*Nx;
M = sparse(iix(1:kx), jjx(1:kx), sx(1:kx), Nx+2, Nx+2);
