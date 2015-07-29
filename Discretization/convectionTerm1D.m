function M = convectionTerm1D(u)
% This function uses the central difference scheme to discretize a 1D
% convection term in the form \grad (u \phi) where u is a face vactor
% It also returns the x and y parts of the matrix of coefficient.
%
% SYNOPSIS:
%  M = convectionTerm1D(u)
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
G = [1:Nx+2];
DXe = u.domain.cellsize.x(3:end);
DXw = u.domain.cellsize.x(1:end-2);
DXp = u.domain.cellsize.x(2:end-1);

% define the vectors to stores the sparse matrix data
iix = zeros(3*(Nx+2),1);
jjx = zeros(3*(Nx+2),1);
sx = zeros(3*(Nx+2),1);

% reassign the east, west for code readability
ue = u.xvalue(2:Nx+1)./(DXp+DXe);
uw = u.xvalue(1:Nx)./(DXp+DXw);

% calculate the coefficients for the internal cells
AE = reshape(ue,Nx,1);
AW = reshape(-uw,Nx,1);
APx = reshape((ue.*DXe-uw.*DXw)./DXp,Nx,1);

% build the sparse matrix based on the numbering system
rowx_index = reshape(G(2:Nx+1),Nx,1); % main diagonal x
iix(1:3*Nx) = repmat(rowx_index,3,1);
jjx(1:3*Nx) = [reshape(G(1:Nx),Nx,1); ...
		reshape(G(2:Nx+1),Nx,1); reshape(G(3:Nx+2),Nx,1)];
sx(1:3*Nx) = [AW; APx; AE];

% build the sparse matrix
kx = 3*Nx;
M = sparse(iix(1:kx), jjx(1:kx), sx(1:kx), Nx+2, Nx+2);
