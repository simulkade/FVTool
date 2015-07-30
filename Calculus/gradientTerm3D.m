function faceGrad = gradientTerm3D(phi)
% this function calculates the gradient of a variable in x and y direction
% it checks for the availability of the ghost variables and use them, otherwise
% estimate them, assuming a zero gradient on the boundaries
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

% check the size of the variable and the mesh dimension
Nx = phi.domain.dims(1);
Ny = phi.domain.dims(2);
Nz = phi.domain.dims(3);
DX = repmat(phi.domain.cellsize.x, 1, Ny, Nz);
DY = repmat(phi.domain.cellsize.y', Nx, 1, Nz);
DZ = zeros(1,1,Nz+2);
DZ(1,1,:) = phi.domain.cellsize.z;
DZ=repmat(DZ, Nx, Ny, 1);
dx = 0.5*(DX(1:end-1,:,:)+DX(2:end,:,:));
dy = 0.5*(DY(:,1:end-1,:)+DY(:,2:end,:));
dz = 0.5*(DZ(:,:,1:end-1)+DZ(:,:,2:end));


% in this case, ghost cells have values
xvalue = (phi.value(2:Nx+2,2:Ny+1,2:Nz+1)-phi.value(1:Nx+1,2:Ny+1,2:Nz+1))./dx;
yvalue = (phi.value(2:Nx+1,2:Ny+2,2:Nz+1)-phi.value(2:Nx+1,1:Ny+1,2:Nz+1))./dy;
zvalue = (phi.value(2:Nx+1,2:Ny+1,2:Nz+2)-phi.value(2:Nx+1,2:Ny+1,1:Nz+1))./dz;


faceGrad=FaceVariable(phi.domain, xvalue, yvalue, zvalue);

