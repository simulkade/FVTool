function M = convectionTermSpherical1D(u)
% This function uses the central difference scheme to discretize a 1D
% convection term in the form \grad (u \phi) where u is a face vactor
% It is for a cylindrical coordinate in the r direction
%
% SYNOPSIS:
%   M = convectionTermCylindrical1D(u)
%
% PARAMETERS:
%	u   - FaceVariable  
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
G = 1:Nr+2;
DXe = u.domain.cellsize.x(3:end);
DXw = u.domain.cellsize.x(1:end-2);
DXp = u.domain.cellsize.x(2:end-1);
% rp = u.domain.cellcenters.x;
rf = u.domain.facecenters.x;

% define the vectors to stores the sparse matrix data
iix = zeros(3*(Nr+2),1);
jjx = zeros(3*(Nr+2),1);
sx = zeros(3*(Nr+2),1);

% reassign the east, west, north, and south velocity vectors for the
% code readability
ue = u.xvalue(2:Nr+1).*DXp./(DXp+DXe).*rf(2:Nr+1).^2./(1/3*(rf(2:Nr+1).^3-rf(1:Nr).^3));
uw = u.xvalue(1:Nr).*DXp./(DXp+DXw).*rf(1:Nr).^2./(1/3*(rf(2:Nr+1).^3-rf(1:Nr).^3));

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
