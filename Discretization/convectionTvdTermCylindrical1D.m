function [M, RHS] = convectionTvdTermCylindrical1D(MeshStructure, u, phi, FL)
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
G = MeshStructure.numbering;
Nr = MeshStructure.numberofcells;
dx = MeshStructure.cellsize;
r = MeshStructure.cellcenters.x';
rf = MeshStructure.facecenters.x';
RHS = zeros(Nr+2, 1);
psi_p = zeros(Nr+1,1);
psi_m = zeros(Nr+1,1);

% define the vectors to store the sparse matrix data
iix = zeros(3*(Nr+2),1);
jjx = zeros(3*(Nr+2),1);
sx = zeros(3*(Nr+2),1);

% extract the velocity data 
% note: size(ux) = [1:m+1, 1:n] and size(uy) = [1:m, 1:n+1]
ux = u.xvalue;

% calculate the upstream to downstream gradient ratios for u>0 (+ ratio)
% P is 3:Nr+2
% W is 2:Nr+1
% WW is 1:Nr
dphi_p = phi(2:Nr+2)-phi(1:Nr+1);
rp = dphi_p(1:end-1)./fsign(dphi_p(2:end));
psi_p(2:Nr+1) = 0.5*FL(rp).*(phi(3:Nr+2)-phi(2:Nr+1));
psi_p(1) = 0; % left boundary will be handled explicitly

% calculate the upstream to downstream gradient ratios for u<0 (- ratio)
% P is 3:Nr+2
% W is 2:Nr+1
% WW is 1:Nr
rm = dphi_p(2:end)./fsign(dphi_p(1:end-1));
psi_m(1:Nr) = 0.5*FL(rm).*(phi(1:Nr)-phi(2:Nr+1));
psi_m(Nr+1) = 0; % right boundary will be handled explicitly

% reassign the east, west, north, and south velocity vectors for the 
% code readability
ue = ux(2:Nr+1);		uw = ux(1:Nr);
re = rf(2:Nr+1);     rw = rf(1:Nr);

% find the velocity direction for the upwind scheme
ue_min = min(ue,0);	ue_max = max(ue,0);
uw_min = min(uw,0);	uw_max = max(uw,0);

% calculate the TVD correction term
RHS(2:Nr+1) = -(1./(dx*r)).*(re.*(ue_max.*psi_p(2:Nr+1)+ue_min.*psi_m(2:Nr+1))- ...
              rw.*(uw_max.*psi_p(1:Nr)+uw_min.*psi_m(1:Nr)));

% calculate the coefficients for the internal cells
AE = reshape(re.*ue_min./(r*dx),Nr,1);
AW = reshape(-rw.*uw_max./(r*dx),Nr,1);
APx = reshape((re.*ue_max-rw.*uw_min)./(r*dx),Nr,1);

% correct for the cells next to the boundary
% Left boundary:
APx(1) = APx(1)-rw(1)*uw_max(1)/(2*r(1)*dx);   AW(1) = AW(1)/2;
% Right boundary:
AE(end) = AE(end)/2;    APx(end) = APx(end) + re(end)*ue_min(end)/(2*r(end)*dx);

% build the sparse matrix based on the numbering system
rowx_index = reshape(G(2:Nr+1),Nr,1); % main diagonal x
iix(1:3*Nr) = repmat(rowx_index,3,1);
jjx(1:3*Nr) = [reshape(G(1:Nr),Nr,1); reshape(G(2:Nr+1),Nr,1); reshape(G(3:Nr+2),Nr,1)];
sx(1:3*Nr) = [AW; APx; AE];

% build the sparse matrix
kx = 3*Nr;
M = sparse(iix(1:kx), jjx(1:kx), sx(1:kx), Nr+2, Nr+2);
end

function phi_out = fsign(phi_in)
% This function checks the value of phi_in and assigns an eps value to the
% elements that are less than or equal to zero, while keeping the signs of
% the nonzero elements
    phi_out = (abs(phi_in)>=eps).*phi_in+eps*(phi_in==0)+eps*(abs(phi_in)<eps).*sign(phi_in);
end
