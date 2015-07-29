function [M, RHS, Mx, My, RHSx, RHSy] = ...
    convectionTvdTerm2D(u, phi, FL)
% This function uses the upwind scheme to discretize a 1D
% convection term in the form \grad (u \phi) where u is a face vactor
% It also returns the x and y parts of the matrix of coefficient.
%
% SYNOPSIS:
%   [M, RHS, Mx, My, RHSx, RHSy] = ...
%    convectionTvdTerm2D(u, phi, FL)
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

%{
Copyright (c) 2012, 2013, Ali Akbar Eftekhari
All rights reserved.

Redistribution and use in source and binary forms, with or
without modification, are permitted provided that the following
conditions are met:

    *   Redistributions of source code must retain the above copyright notice,
        this list of conditions and the following disclaimer.
    *   Redistributions in binary form must reproduce the above
        copyright notice, this list of conditions and the following
        disclaimer in the documentation and/or other materials provided
        with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR
CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA,
OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
%}

% extract data from the mesh structure
Nx = u.domain.dims(1);
Ny = u.domain.dims(2);
G=reshape(1:(Nx+2)*(Ny+2), Nx+2, Ny+2);
DXp = repmat(u.domain.cellsize.x(2:end-1), 1, Ny);
DYp = repmat(u.domain.cellsize.y(2:end-1)', Nx, 1);
dx=repmat(0.5*(u.domain.cellsize.x(1:end-1)+u.domain.cellsize.x(2:end)), 1, Ny);
dy=repmat(0.5*(u.domain.cellsize.y(1:end-1)+u.domain.cellsize.y(2:end))', Nx, 1);
psiX_p = zeros(Nx+1,Ny);
psiX_m = zeros(Nx+1,Ny);
psiY_p = zeros(Nx,Ny+1);
psiY_m = zeros(Nx,Ny+1);

% define the vectors to stores the sparse matrix data
iix = zeros(3*(Nx+2)*(Ny+2),1);	iiy = zeros(3*(Nx+2)*(Ny+2),1);
jjx = zeros(3*(Nx+2)*(Ny+2),1);	jjy = zeros(3*(Nx+2)*(Ny+2),1);
sx = zeros(3*(Nx+2)*(Ny+2),1);	sy = zeros(3*(Nx+2)*(Ny+2),1);
mnx = Nx*Ny;	mny = Nx*Ny;

% extract the velocity data
% note: size(ux) = [1:m+1, 1:n] and size(uy) = [1:m, 1:n+1]
ux = u.xvalue;
uy = u.yvalue;

% calculate the upstream to downstream gradient ratios for u>0 (+ ratio)
% x direction
dphiX_p = (phi.value(2:Nx+2, 2:Ny+1)-phi.value(1:Nx+1, 2:Ny+1))./dx;
rX_p = dphiX_p(1:end-1,:)./fsign(dphiX_p(2:end,:));
psiX_p(2:Nx+1,:) = 0.5*FL(rX_p).*(phi(3:Nx+2,2:Ny+1)-phi(2:Nx+1, 2:Ny+1));
psiX_p(1, :) = 0; % left boundary will be handled in the main matrix
% y direction
dphiY_p = (phi(2:Nx+1, 2:Ny+2)-phi(2:Nx+1, 1:Ny+1))./dy;
rY_p = dphiY_p(:,1:end-1)./fsign(dphiY_p(:,2:end));
psiY_p(:,2:Ny+1) = 0.5*FL(rY_p).*(phi(2:Nx+1,3:Ny+2)-phi(2:Nx+1, 2:Ny+1));
psiY_p(:,1) = 0; % Bottom boundary will be handled in the main matrix

% calculate the upstream to downstream gradient ratios for u<0 (- ratio)
% x direction
rX_m = dphiX_p(2:end,:)./fsign(dphiX_p(1:end-1,:));
psiX_m(1:Nx,:) = 0.5*FL(rX_m).*(phi(1:Nx, 2:Ny+1)-phi(2:Nx+1, 2:Ny+1));
psiX_m(Nx+1,:) = 0; % right boundary
% y direction
rY_m = dphiY_p(:,2:end)./fsign(dphiY_p(:,1:end-1));
psiY_m(:,1:Ny) = 0.5*FL(rY_m).*(phi(2:Nx+1, 1:Ny)-phi(2:Nx+1, 2:Ny+1));
psiY_m(:, Ny+1) = 0; % top boundary will be handled in the main matrix

% reassign the east, west, north, and south velocity vectors for the
% code readability
ue = ux(2:Nx+1,:);		uw = ux(1:Nx,:);
vn = uy(:,2:Ny+1);       vs = uy(:,1:Ny);

% find the velocity direction for the upwind scheme
ue_min = min(ue,0);	ue_max = max(ue,0);
uw_min = min(uw,0);	uw_max = max(uw,0);
vn_min = min(vn,0);	vn_max = max(vn,0);
vs_min = min(vs,0);	vs_max = max(vs,0);

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

% calculate the TVD correction term
div_x = -(1./DXp).*((ue_max.*psiX_p(2:Nx+1,:)+ue_min.*psiX_m(2:Nx+1,:))- ...
              (uw_max.*psiX_p(1:Nx,:)+uw_min.*psiX_m(1:Nx,:)));
div_y = -(1./DYp).*((vn_max.*psiY_p(:,2:Ny+1)+vn_min.*psiY_m(:,2:Ny+1))- ...
              (vs_max.*psiY_p(:,1:Ny)+vs_min.*psiY_m(:,1:Ny)));
% define the RHS Vector
RHS = zeros((Nx+2)*(Ny+2),1);
RHSx = zeros((Nx+2)*(Ny+2),1);
RHSy = zeros((Nx+2)*(Ny+2),1);

% assign the values of the RHS vector
RHS(rowx_index) = reshape(div_x+div_y,Nx*Ny,1);
RHSx(rowx_index) = reshape(div_x,Nx*Ny,1);
RHSy(rowy_index) = reshape(div_y,Nx*Ny,1);

% build the sparse matrix
kx = 3*mnx;
ky = 3*mny;
Mx = sparse(iix(1:kx), jjx(1:kx), sx(1:kx), (Nx+2)*(Ny+2), (Nx+2)*(Ny+2));
My = sparse(iiy(1:ky), jjy(1:ky), sy(1:ky), (Nx+2)*(Ny+2), (Nx+2)*(Ny+2));
M = Mx + My;

end

function phi_out = fsign(phi_in)
% This function checks the value of phi_in and assigns an eps value to the
% elements that are less than or equal to zero, while keeping the signs of
% the nonzero elements
    phi_out = (abs(phi_in)>=eps).*phi_in+eps*(phi_in==0)+eps*(abs(phi_in)<eps).*sign(phi_in);
end
