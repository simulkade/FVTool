function RHSout = excludeGhostRHS(MS, RHS)
% this function cuts out the RHS values related to ghost cells and returns
% an RHS, which is only for internal cells
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
d = MS.dimension;

if (d ==1) || (d==1.5)
	Nx=MS.dims(1);
	G = 1:Nx+2;
	RHSout = RHS(reshape(G(2:end-1),Nx,1));
elseif (d == 2) || (d == 2.5) || (d==2.8)
	Nxy = MS.dims;
	Nx = Nxy(1); Ny = Nxy(2);
	G=reshape(1:(Nx+2)*(Ny+2), Nx+2, Ny+2);
	RHSout = RHS(reshape(G(2:end-1,2:end-1),Nx*Ny,1));
elseif (d == 3) || (d==3.2)
	Nxyz = MS.dims;
	Nx = Nxyz(1); Ny = Nxyz(2); Nz = Nxyz(3);
	G=reshape(1:(Nx+2)*(Ny+2)*(Nz+2), Nx+2, Ny+2, Nz+2);
  RHSout = RHS(reshape(G(2:end-1,2:end-1,2:end-1),Nx*Ny*Nz,1));
end
