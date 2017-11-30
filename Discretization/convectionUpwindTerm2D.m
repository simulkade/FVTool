function [M, Mx, My] = convectionUpwindTerm2D(u, varargin)
% This function uses the upwind scheme to discretize a 2D
% convection term in the form \grad (u \phi) where u is a face vactor
% It also returns the x and y parts of the matrix of coefficient.
%
% SYNOPSIS:
%   (M, Mx, My) = convectionUpwindTerm2D(u)
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
Ny = u.domain.dims(2);
G=reshape(1:(Nx+2)*(Ny+2), Nx+2, Ny+2);
DXp = repmat(u.domain.cellsize.x(2:end-1), 1, Ny);
DYp = repmat(u.domain.cellsize.y(2:end-1)', Nx, 1);

if nargin>1
    u_upwind = varargin{1};
else
    u_upwind = u;
end

% define the vectors to stores the sparse matrix data
iix = zeros(3*(Nx+2)*(Ny+2),1);	iiy = zeros(3*(Nx+2)*(Ny+2),1);
jjx = zeros(3*(Nx+2)*(Ny+2),1);	jjy = zeros(3*(Nx+2)*(Ny+2),1);
sx = zeros(3*(Nx+2)*(Ny+2),1);	sy = zeros(3*(Nx+2)*(Ny+2),1);
mnx = Nx*Ny;	mny = Nx*Ny;

ue_min = u.xvalue(2:Nx+1,:);
ue_max = u.xvalue(2:Nx+1,:);
uw_min = u.xvalue(1:Nx,:);
uw_max = u.xvalue(1:Nx,:);
vn_min = u.yvalue(:,2:Ny+1);
vn_max = u.yvalue(:,2:Ny+1);
vs_min = u.yvalue(:,1:Ny);
vs_max = u.yvalue(:,1:Ny);

ue_min(u_upwind.xvalue(2:Nx+1,:)>0.0) = 0.0;
ue_max(u_upwind.xvalue(2:Nx+1,:)<0.0) = 0.0;
uw_min(u_upwind.xvalue(1:Nx,:)>0.0) = 0.0;
uw_max(u_upwind.xvalue(1:Nx,:)<0.0) = 0.0;
vn_min(u_upwind.yvalue(:,2:Ny+1)>0.0) = 0.0;
vn_max(u_upwind.yvalue(:,2:Ny+1)<0.0) = 0.0;
vs_min(u_upwind.yvalue(:,1:Ny)>0.0) = 0.0;
vs_max(u_upwind.yvalue(:,1:Ny)<0.0) = 0.0;

% calculate the coefficients for the internal cells
AE = ue_min./DXp;
AW = -uw_max./DXp;
AN = vn_min./DYp;
AS = -vs_max./DYp;
APx = (ue_max-uw_min)./DXp;
APy = (vn_max-vs_min)./DYp;

% Also correct for the boundary cells (not the ghost cells)
% Left boundary:
APx(1,:) = APx(1,:)-uw_max(1,:)/(2*DXp(1));   AW(1,:) = AW(1,:)/2;
% Right boundary:
AE(end,:) = AE(end,:)/2;    APx(end,:) = APx(end,:) + ue_min(end,:)/(2*DXp(end));
% Bottom boundary:
APy(:,1) = APy(:,1)-vs_max(:,1)/(2*DYp(1));   AS(:,1) = AS(:,1)/2;
% Top boundary:
AN(:,end) = AN(:,end)/2;    APy(:,end) = APy(:,end) + vn_min(:,end)/(2*DYp(end));

AE = reshape(AE,mnx,1);
AW = reshape(AW,mnx,1);
AN = reshape(AN,mny,1);
AS = reshape(AS,mny,1);
APx = reshape(APx,mnx,1);
APy = reshape(APy,mny,1);

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
