function M = convectionTermCylindrical1D(u)
% This function uses the central difference scheme to discretize a 1D
% convection term in the form \grad (u \phi) where u is a face vactor
% It is for a cylindrical coordinate in the r direction
%
% SYNOPSIS:
%   M = convectionTermCylindrical1D(u)
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
Nr = u.domain.dims(1);
G = 1:Nr+2;
DXe = u.domain.cellsize.x(3:end);
DXw = u.domain.cellsize.x(1:end-2);
DXp = u.domain.cellsize.x(2:end-1);
rp = u.domain.cellcenters.x;
rf = u.domain.facecenters.x;

% define the vectors to stores the sparse matrix data
iix = zeros(3*(Nr+2),1);
jjx = zeros(3*(Nr+2),1);
sx = zeros(3*(Nr+2),1);

% extract the velocity data
% note: size(ux) = [1:m+1, 1:n] and size(uy) = [1:m, 1:n+1]
ux = u.xvalue;

% reassign the east, west, north, and south velocity vectors for the
% code readability
ue = rf(2:Nr+1).*u.xvalue(2:Nr+1)./(rp.*(DXp+DXe));
uw = rf(1:Nr).*u.xvalue(1:Nr)./(rp.*(DXp+DXw));

% calculate the coefficients for the internal cells
AE = reshape(ue,Nr,1);
AW = reshape(-uw,Nr,1);
APx = reshape((ue.*DXe-uw.*DXw)./DXp,Nr,1);

% build the sparse matrix based on the numbering system
rowx_index = reshape(G(2:Nr+1),Nr,1); % main diagonal x
iix(1:3*Nr) = repmat(rowx_index,3,1);
jjx(1:3*Nr) = [reshape(G(1:Nr),Nr,1); ...
		reshape(G(2:Nr+1),Nr,1); reshape(G(3:Nr+2),Nr,1)];
sx(1:3*Nr) = [AW; APx; AE];

% build the sparse matrix
kx = 3*Nr;
M = sparse(iix(1:kx), jjx(1:kx), sx(1:kx), Nr+2, Nr+2);
