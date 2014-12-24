function faceGrad = gradientTermCylindrical3D(MeshStructure, phi)
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
Nxyz = MeshStructure.numberofcells;
Nx = Nxyz(1); Ntetta = Nxyz(2); Nz = Nxyz(3);
Nxyz_phi = size(phi);
Nx_phi = Nxyz_phi(1);
d = MeshStructure.cellsize;
dx = d(1); dtetta = d(2); dz = d(3);

if Nx_phi == Nx
    % define the average face variabe
    rp = repmat(MeshStructure.cellcenters.x', 1, Ntetta-1, Nz);
    faceGrad.xvalue = zeros(Nx+1,Ntetta,Nz);
    faceGrad.yvalue = zeros(Nx,Ntetta+1,Nz);
    faceGrad.zvalue = zeros(Nx,Ntetta,Nz+1);
    faceGrad.xvalue(2:Nx,:,:) = (phi(2:Nx,:,:)-phi(1:Nx-1,:,:))/dx;
    faceGrad.yvalue(:,2:Ntetta,:) = (phi(:,2:Ntetta,:)-phi(:,1:Ntetta-1,:))./(dtetta*rp);
    faceGrad.zvalue(:,:,2:Nz) = (phi(:,:,2:Nz)-phi(:,:,1:Nz-1))/dz;    
else
    % in this case, ghost cells have values
    rp = repmat(MeshStructure.cellcenters.x', 1, Ntetta+1, Nz);
    faceGrad.xvalue = (phi(2:Nx+2,2:Ntetta+1,2:Nz+1)-phi(1:Nx+1,2:Ntetta+1,2:Nz+1))/dx;
    faceGrad.yvalue = (phi(2:Nx+1,2:Ntetta+2,2:Nz+1)-phi(2:Nx+1,1:Ntetta+1,2:Nz+1))./(dtetta*rp);
    faceGrad.zvalue = (phi(2:Nx+1,2:Ntetta+1,2:Nz+2)-phi(2:Nx+1,2:Ntetta+1,1:Nz+1))/dz;
end


