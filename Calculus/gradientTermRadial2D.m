function faceGrad = gradientTermRadial2D(MeshStructure, phi)
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
Nxy = MeshStructure.numberofcells;
Nx = Nxy(1); Ntetta = Nxy(2);
Nxy_phi = size(phi);
Nx_phi = Nxy_phi(1);
dxdy = MeshStructure.cellsize;
dx = dxdy(1); dtetta = dxdy(2);
if Nx_phi == Nx
    rp = repmat(MeshStructure.cellcenters.x', 1, Ntetta-1);
    % define the average face variabe
    faceGrad.xvalue = [zeros(1,Ntetta); (phi(2:Nx,:)-phi(1:Nx-1,:))/dx; zeros(1,Ntetta)];
    faceGrad.yvalue = [zeros(Nx,1), (phi(:,2:Ntetta)-phi(:,1:Ntetta-1))./(dtetta*rp), zeros(Nx,1)];
else
    rp = repmat(MeshStructure.cellcenters.x', 1, Ntetta+1);
    % in this case, ghost cells have values
    faceGrad.xvalue = (phi(2:Nx+2,2:Ntetta+1)-phi(1:Nx+1,2:Ntetta+1))/dx;
    faceGrad.yvalue = (phi(2:Nx+1,2:Ntetta+2)-phi(2:Nx+1,1:Ntetta+1))./(dtetta*rp);
end


