function [M, Mx, My, Mz] = convectionTermCylindrical3D(u)
% This function uses the central difference scheme to discretize a 2D
% convection term in the form \grad (u \phi) where u is a face vactor
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
Nr = u.domain.dims(1);
Ntetta = u.domain.dims(2);
Nz = u.domain.dims(3);
G=reshape((1:(Nr+2)*(Ntetta+2)*(Nz+2)), Nr+2, Ntetta+2, Nz+2);
DRe = repmat(u.domain.cellsize.x(3:end), 1, Ntetta, Nz);
DRw = repmat(u.domain.cellsize.x(1:end-2), 1, Ntetta, Nz);
DRp = repmat(u.domain.cellsize.x(2:end-1), 1, Ntetta, Nz);
DTHETAn = repmat(u.domain.cellsize.y(3:end)', Nr, 1, Nz);
DTHETAs = repmat(u.domain.cellsize.y(1:end-2)', Nr, 1, Nz);
DTHETAp = repmat(u.domain.cellsize.y(2:end-1)', Nr, 1, Nz);
DZ = zeros(1,1,Nz+2);
DZ(1,1,:) = u.domain.cellsize.z;
DZf=repmat(DZ(1,1,3:end), Nr, Ntetta, 1);
DZb=repmat(DZ(1,1,1:end-2), Nr, Ntetta, 1);
DZp=repmat(DZ(1,1,2:end-1), Nr, Ntetta, 1);
rp = repmat(u.domain.cellcenters.x, 1, Ntetta, Nz);
rf = repmat(u.domain.facecenters.x, 1, Ntetta, Nz);

% define the vectors to stores the sparse matrix data
iix = zeros(3*(Nr+2)*(Ntetta+2)*(Nz+2),1);
jjx = zeros(3*(Nr+2)*(Ntetta+2)*(Nz+2),1);
sx = zeros(3*(Nr+2)*(Ntetta+2)*(Nz+2),1);
iiy = zeros(3*(Nr+2)*(Ntetta+2)*(Nz+2),1);
jjy = zeros(3*(Nr+2)*(Ntetta+2)*(Nz+2),1);
sy = zeros(3*(Nr+2)*(Ntetta+2)*(Nz+2),1);
iiz = zeros(3*(Nr+2)*(Ntetta+2)*(Nz+2),1);
jjz = zeros(3*(Nr+2)*(Ntetta+2)*(Nz+2),1);
sz = zeros(3*(Nr+2)*(Ntetta+2)*(Nz+2),1);
mnx = Nr*Ntetta*Nz;	mny = Nr*Ntetta*Nz;   mnz = Nr*Ntetta*Nz;

% reassign the east, west, north, and south velocity vectors for the
% code readability
ue = rf(2:Nr+1,:,:).*u.xvalue(2:Nr+1,:,:)./(rp.*(DRp+DRe));
uw = rf(1:Nr,:,:).*u.xvalue(1:Nr,:,:)./(rp.*(DRp+DRw));
vn = u.yvalue(:,2:Ntetta+1,:)./(rp.*(DTHETAp+DTHETAn));
vs = u.yvalue(:,1:Ntetta,:)./(rp.*(DTHETAp+DTHETAs));
wf = u.zvalue(:,:,2:Nz+1)./(DZp+DZf);
wb = u.zvalue(:,:,1:Nz)./(DZp+DZb);

% calculate the coefficients for the internal cells
AE = reshape(ue,mnx,1);
AW = reshape(-uw,mnx,1);
AN = reshape(vn,mny,1);
AS = reshape(-vs,mny,1);
AF = reshape(wf,mnz,1);
AB = reshape(-wb,mnz,1);
APx = reshape((DRe.*ue-DRw.*uw)./DRp,mnx,1);
APy = reshape((DTHETAn.*vn-DTHETAs.*vs)./DTHETAp,mny,1);
APz = reshape((DZf.*wf-DZb.*wb)./DZp,mnz,1);

% build the sparse matrix based on the numbering system
rowx_index = reshape(G(2:Nr+1,2:Ntetta+1,2:Nz+1),mnx,1); % main diagonal x
iix(1:3*mnx) = repmat(rowx_index,3,1);
rowy_index = reshape(G(2:Nr+1,2:Ntetta+1,2:Nz+1),mny,1); % main diagonal y
iiy(1:3*mny) = repmat(rowy_index,3,1);
rowz_index = reshape(G(2:Nr+1,2:Ntetta+1,2:Nz+1),mnz,1); % main diagonal z
iiz(1:3*mnz) = repmat(rowz_index,3,1);
jjx(1:3*mnx) = [reshape(G(1:Nr,2:Ntetta+1,2:Nz+1),mnx,1); reshape(G(2:Nr+1,2:Ntetta+1,2:Nz+1),mnx,1); reshape(G(3:Nr+2,2:Ntetta+1,2:Nz+1),mnx,1)];
jjy(1:3*mny) = [reshape(G(2:Nr+1,1:Ntetta,2:Nz+1),mny,1); reshape(G(2:Nr+1,2:Ntetta+1,2:Nz+1),mny,1); reshape(G(2:Nr+1,3:Ntetta+2,2:Nz+1),mny,1)];
jjz(1:3*mnz) = [reshape(G(2:Nr+1,2:Ntetta+1,1:Nz),mnz,1); reshape(G(2:Nr+1,2:Ntetta+1,2:Nz+1),mnz,1); reshape(G(2:Nr+1,2:Ntetta+1,3:Nz+2),mnz,1)];
sx(1:3*mnx) = [AW; APx; AE];
sy(1:3*mny) = [AS; APy; AN];
sz(1:3*mnz) = [AB; APz; AF];

% build the sparse matrix
kx = 3*mnx;
ky = 3*mny;
kz = 3*mnz;
Mx = sparse(iix(1:kx), jjx(1:kx), sx(1:kx), (Nr+2)*(Ntetta+2)*(Nz+2), (Nr+2)*(Ntetta+2)*(Nz+2));
My = sparse(iiy(1:ky), jjy(1:ky), sy(1:ky), (Nr+2)*(Ntetta+2)*(Nz+2), (Nr+2)*(Ntetta+2)*(Nz+2));
Mz = sparse(iiz(1:kz), jjz(1:kz), sz(1:kz), (Nr+2)*(Ntetta+2)*(Nz+2), (Nr+2)*(Ntetta+2)*(Nz+2));
M = Mx + My + Mz;
