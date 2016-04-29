function M = convectionUpwindTerm1D(u)
% This function uses the upwind scheme to discretize a 1D
% convection term in the form \grad (u \phi) where u is a face vactor
% It also returns the x and y parts of the matrix of coefficient.
%
% SYNOPSIS:
%   M = convectionUpwindTerm1D(u)
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
Nx = u.domain.dims(1);
G = [1:Nx+2];
DXp = u.domain.cellsize.x(2:end-1);

% define the vectors to store the sparse matrix data
iix = zeros(3*(Nx+2),1);
jjx = zeros(3*(Nx+2),1);
sx = zeros(3*(Nx+2),1);

% extract the velocity data
% note: size(ux) = [1:m+1, 1:n] and size(uy) = [1:m, 1:n+1]
ux = u.xvalue;

% reassign the east, west, north, and south velocity vectors for the
% code readability
ue = ux(2:Nx+1);		uw = ux(1:Nx);

% find the velocity direction for the upwind scheme
ue_min = min(ue,0);	ue_max = max(ue,0);
uw_min = min(uw,0);	uw_max = max(uw,0);

% calculate the coefficients for the internal cells
AE = reshape(ue_min./DXp,Nx,1);
AW = reshape(-uw_max./DXp,Nx,1);
APx = reshape((ue_max-uw_min)./DXp,Nx,1);

% correct for the cells next to the boundary
% Left boundary:
APx(1) = APx(1)-uw_max(1)/(2*DXp(1));   AW(1) = AW(1)/2;
% Right boundary:
AE(end) = AE(end)/2;    APx(end) = APx(end) + ue_min(end)/(2*DXp(end));

% build the sparse matrix based on the numbering system
rowx_index = reshape(G(2:Nx+1),Nx,1); % main diagonal x
iix(1:3*Nx) = repmat(rowx_index,3,1);
jjx(1:3*Nx) = [reshape(G(1:Nx),Nx,1); reshape(G(2:Nx+1),Nx,1); reshape(G(3:Nx+2),Nx,1)];
sx(1:3*Nx) = [AW; APx; AE];

% build the sparse matrix
kx = 3*Nx;
M = sparse(iix(1:kx), jjx(1:kx), sx(1:kx), Nx+2, Nx+2);
