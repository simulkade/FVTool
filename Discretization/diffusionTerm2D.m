function [M, Mx, My] = diffusionTerm2D(D)
% This function uses the central difference scheme to discretize a 2D
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
Ny = D.domain.dims(2);
G=reshape(1:(Nx+2)*(Ny+2), Nx+2, Ny+2);

DX = repmat(D.domain.cellsize.x, 1, Ny);
DY = repmat(D.domain.cellsize.y', Nx, 1);
dx = 0.5*(DX(1:end-1,:)+DX(2:end,:));
dy = 0.5*(DY(:,1:end-1)+DY(:,2:end));

% define the vectors to store the sparse matrix data
iix = zeros(3*(Nx+2)*(Ny+2),1);	iiy = zeros(3*(Nx+2)*(Ny+2),1);
jjx = zeros(3*(Nx+2)*(Ny+2),1);	jjy = zeros(3*(Nx+2)*(Ny+2),1);
sx = zeros(3*(Nx+2)*(Ny+2),1);	sy = zeros(3*(Nx+2)*(Ny+2),1);
mnx = Nx*Ny;	mny = Nx*Ny;

% reassign the east, west, north, and south velocity vectors for the
% code readability
De = D.xvalue(2:Nx+1,:)./(dx(2:Nx+1,:).*DX(2:Nx+1,:));
Dw = D.xvalue(1:Nx,:)./(dx(1:Nx,:).*DX(2:Nx+1,:));
Dn = D.yvalue(:,2:Ny+1)./(dy(:,2:Ny+1).*DY(:,2:Ny+1));
Ds = D.yvalue(:,1:Ny)./(dy(:,1:Ny).*DY(:,2:Ny+1));

% calculate the coefficients for the internal cells
AE = reshape(De,mnx,1);
AW = reshape(Dw,mnx,1);
AN = reshape(Dn,mny,1);
AS = reshape(Ds,mny,1);
APx = reshape(-(De+Dw),mnx,1);
APy = reshape(-(Dn+Ds),mny,1);

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
