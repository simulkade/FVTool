function [M, RHS, Mx, My, RHSx, RHSy] = ...
    convectionTvdTermCylindrical2D(MeshStructure, u, phi, FL)
% This function uses the upwind scheme to discretize a 1D
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
AND ANz EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, 
THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR 
PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR 
CONTRIBUTORS BE LIABLE FOR ANz DIRECT, INDIRECT, INCIDENTAL, SPECIAL, 
EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, 
PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, 
OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANz THEORY OF 
LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING 
NEGLIGENCE OR OTHERWISE) ARISING IN ANz WAY OUT OF THE USE OF THIS 
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
%}

% extract data from the mesh structure
G = MeshStructure.numbering;
Nrz = MeshStructure.numberofcells;
Nr = Nrz(1); Nz = Nrz(2);
dxdy = MeshStructure.cellsize;
dx = dxdy(1); dy = dxdy(2);
r = repmat(MeshStructure.cellcenters.x', 1, Nz);
rf = repmat(MeshStructure.facecenters.x', 1, Nz);
psiX_p = zeros(Nr+1,Nz);
psiX_m = zeros(Nr+1,Nz);
psiY_p = zeros(Nr,Nz+1);
psiY_m = zeros(Nr,Nz+1);

% define the vectors to stores the sparse matrix data
iix = zeros(3*(Nr+2)*(Nz+2),1);	iiy = zeros(3*(Nr+2)*(Nz+2),1);
jjx = zeros(3*(Nr+2)*(Nz+2),1);	jjy = zeros(3*(Nr+2)*(Nz+2),1);
sx = zeros(3*(Nr+2)*(Nz+2),1);	sy = zeros(3*(Nr+2)*(Nz+2),1);
mNr = Nr*Nz;	mNz = Nr*Nz;

% extract the velocity data 
% note: size(ux) = [1:m+1, 1:n] and size(uy) = [1:m, 1:n+1]
ux = u.xvalue;
uy = u.yvalue;

% calculate the upstream to downstream gradient ratios for u>0 (+ ratio)
% x direction
dphiX_p = phi(2:Nr+2, 2:Nz+1)-phi(1:Nr+1, 2:Nz+1);
rX_p = dphiX_p(1:end-1,:)./fsign(dphiX_p(2:end,:));
psiX_p(2:Nr+1,:) = 0.5*FL(rX_p).*(phi(3:Nr+2,2:Nz+1)-phi(2:Nr+1, 2:Nz+1));
psiX_p(1, :) = 0; % left boundary will be handled in the main matrix
% y direction
dphiY_p = phi(2:Nr+1, 2:Nz+2)-phi(2:Nr+1, 1:Nz+1);
rY_p = dphiY_p(:,1:end-1)./fsign(dphiY_p(:,2:end));
psiY_p(:,2:Nz+1) = 0.5*FL(rY_p).*(phi(2:Nr+1,3:Nz+2)-phi(2:Nr+1, 2:Nz+1));
psiY_p(:,1) = 0; % Bottom boundary will be handled in the main matrix

% calculate the upstream to downstream gradient ratios for u<0 (- ratio)
% x direction
rX_m = dphiX_p(2:end,:)./fsign(dphiX_p(1:end-1,:));
psiX_m(1:Nr,:) = 0.5*FL(rX_m).*(phi(1:Nr, 2:Nz+1)-phi(2:Nr+1, 2:Nz+1));
psiX_m(Nr+1,:) = 0; % right boundary
% y direction
rY_m = dphiY_p(:,2:end)./fsign(dphiY_p(:,1:end-1));
psiY_m(:,1:Nz) = 0.5*FL(rY_m).*(phi(2:Nr+1, 1:Nz)-phi(2:Nr+1, 2:Nz+1));
psiY_m(:, Nz+1) = 0; % top boundary will be handled in the main matrix

% reassign the east, west, north, and south velocity vectors for the 
% code readability
ue = ux(2:Nr+1,:);		uw = ux(1:Nr,:);
vn = uy(:,2:Nz+1);       vs = uy(:,1:Nz);
re = rf(2:Nr+1,:);         rw = rf(1:Nr,:);

% find the velocity direction for the upwind scheme
ue_min = min(ue,0);	ue_max = max(ue,0);
uw_min = min(uw,0);	uw_max = max(uw,0);
vn_min = min(vn,0);	vn_max = max(vn,0);
vs_min = min(vs,0);	vs_max = max(vs,0);

% calculate the coefficients for the internal cells
AE = re.*ue_min./(r*dx);
AW = -rw.*uw_max./(r*dx);
AN = vn_min/dy;
AS = -vs_max/dy;
APx = (re.*ue_max-rw.*uw_min)./(r*dx);
APy = (vn_max-vs_min)/dy;

% Also correct for the boundary cells (not the ghost cells)
% Left boundary:
APx(1,:) = APx(1,:)-rw(1,:).*uw_max(1,:)./(2*r(1,:)*dx);   AW(1,:) = AW(1,:)/2;
% Right boundary:
AE(end,:) = AE(end,:)/2;    APx(end,:) = APx(end,:) + re(end,:).*ue_min(end,:)./(2*r(end,:)*dx);
% Bottom boundary:
APy(:,1) = APy(:,1)-vs_max(:,1)/(2*dy);   AS(:,1) = AS(:,1)/2;
% Top boundary:
AN(:,end) = AN(:,end)/2;    APy(:,end) = APy(:,end) + vn_min(:,end)/(2*dy);

AE = reshape(AE,mNr,1);
AW = reshape(AW,mNr,1);
AN = reshape(AN,mNz,1);
AS = reshape(AS,mNz,1);
APx = reshape(APx,mNr,1);
APy = reshape(APy,mNz,1);

% build the sparse matrix based on the numbering system
rowx_index = reshape(G(2:Nr+1,2:Nz+1),mNr,1); % main diagonal x
iix(1:3*mNr) = repmat(rowx_index,3,1);
rowy_index = reshape(G(2:Nr+1,2:Nz+1),mNz,1); % main diagonal y
iiy(1:3*mNz) = repmat(rowy_index,3,1);
jjx(1:3*mNr) = [reshape(G(1:Nr,2:Nz+1),mNr,1); reshape(G(2:Nr+1,2:Nz+1),mNr,1); reshape(G(3:Nr+2,2:Nz+1),mNr,1)];
jjy(1:3*mNz) = [reshape(G(2:Nr+1,1:Nz),mNz,1); reshape(G(2:Nr+1,2:Nz+1),mNz,1); reshape(G(2:Nr+1,3:Nz+2),mNz,1)];
sx(1:3*mNr) = [AW; APx; AE];
sy(1:3*mNz) = [AS; APy; AN];

% calculate the TVD correction term
div_x = -(1./(dx*r)).*(re.*(ue_max.*psiX_p(2:Nr+1,:)+ue_min.*psiX_m(2:Nr+1,:))- ...
              rw.*(uw_max.*psiX_p(1:Nr,:)+uw_min.*psiX_m(1:Nr,:)));
div_y = -(1/dy)*((vn_max.*psiY_p(:,2:Nz+1)+vn_min.*psiY_m(:,2:Nz+1))- ...
              (vs_max.*psiY_p(:,1:Nz)+vs_min.*psiY_m(:,1:Nz)));
% define the RHS Vector
RHS = zeros((Nr+2)*(Nz+2),1);
RHSx = zeros((Nr+2)*(Nz+2),1);
RHSy = zeros((Nr+2)*(Nz+2),1);

% assign the values of the RHS vector
RHS(rowx_index) = reshape(div_x+div_y,Nr*Nz,1);
RHSx(rowx_index) = reshape(div_x,Nr*Nz,1);
RHSy(rowy_index) = reshape(div_y,Nr*Nz,1);

% build the sparse matrix
kx = 3*mNr;
ky = 3*mNz;
Mx = sparse(iix(1:kx), jjx(1:kx), sx(1:kx), (Nr+2)*(Nz+2), (Nr+2)*(Nz+2));
My = sparse(iiy(1:ky), jjy(1:ky), sy(1:ky), (Nr+2)*(Nz+2), (Nr+2)*(Nz+2));
M = Mx + My;

end

function phi_out = fsign(phi_in)
% This function checks the value of phi_in and assigns an eps value to the
% elements that are less than or equal to zero, while keeping the signs of
% the nonzero elements
    phi_out = (abs(phi_in)>=eps).*phi_in+eps*(phi_in==0)+eps*(abs(phi_in)<eps).*sign(phi_in);
end
