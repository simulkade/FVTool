function phiFaceAverage = arithmeticMean3D(MeshStructure, phi)
% This function gets the value of the field variable phi defined
% over the MeshStructure and calculates the arithmetic average on 
% the cell faces, for a uniform mesh.
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
Nx = Nxyz(1); Ny = Nxyz(2); Nz = Nxyz(3);
Nxyz_phi = size(phi);
Nx_phi = Nxyz_phi(1);
if Nx_phi == Nx
	% define the average face variabe
    phiFaceAverage.xvalue = zeros(Nx+1,Ny,Nz);
    phiFaceAverage.xvalue(1,:,:) = phi(1,:,:);
    phiFaceAverage.xvalue(Nx+1,:,:) = phi(Nx,:,:);
    phiFaceAverage.yvalue = zeros(Nx,Ny+1,Nz);
    phiFaceAverage.yvalue(:,1,:) = phi(:,1,:);
    phiFaceAverage.yvalue(:,Ny+1,:) = phi(:,Ny,:);
    phiFaceAverage.zvalue = zeros(Nx,Ny,Nz+1);
    phiFaceAverage.zvalue(:,:,1) = phi(:,:,1);
    phiFaceAverage.zvalue(:,:,Nz+1) = phi(:,:,Nz);
	phiFaceAverage.xvalue(2:Nx,:,:) = (phi(1:Nx-1,:,:)+phi(2:Nx,:,:))/2;
	phiFaceAverage.yvalue(:,2:Ny,:) = (phi(:,1:Ny-1,:)+phi(:,2:Ny,:))/2;
    phiFaceAverage.zvalue(:,:,2:Nz) = (phi(:,:,1:Nz-1)+phi(:,:,2:Nz))/2;
else
	% in this case, ghost cells have values
	phiFaceAverage.xvalue = (phi(1:Nx+1,2:Ny+1,2:Nz+1)+phi(2:Nx+2,2:Ny+1,2:Nz+1))/2;
	phiFaceAverage.yvalue = (phi(2:Nx+1,1:Ny+1,2:Nz+1)+phi(2:Nx+1,2:Ny+2,2:Nz+1))/2;
    phiFaceAverage.zvalue = (phi(2:Nx+1,2:Ny+1,1:Nz+1)+phi(2:Nx+1,2:Ny+1,2:Nz+2))/2;
end
