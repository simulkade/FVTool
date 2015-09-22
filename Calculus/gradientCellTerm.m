function cellGrad = gradientCellTerm(phi)
% this function calculates the gradient of a variable in x direction, in
% the cell center. It needs a cell variable as an input.
%
% SYNOPSIS:
%   cellGrad = gradientCellTerm(phi)
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

phi_face=linearMean(phi);
d = phi.domain.dimension;
if (d ==1) || (d==1.5)
	Nx = phi.domain.dims(1);
	DX = phi.domain.cellsize.x(2:end-1);
	xvalue = (phi_face.xvalue(2:Nx+1)-phi_face.xvalue(1:Nx))./DX;
	yvalue=[];
	zvalue=[];
	cellGrad=FaceVariable(phi.domain, xvalue, yvalue, zvalue);
elseif (d == 2) || (d == 2.5)
	Nx = phi.domain.dims(1);
	Ny = phi.domain.dims(2);
	DX = repmat(phi.domain.cellsize.x(2:end-1), 1, Ny);
	DY = repmat(phi.domain.cellsize.y(2:end-1)', Nx, 1);
	xvalue = (phi_face.xvalue(2:Nx+1,:)-phi_face.xvalue(1:Nx,:))./DX;
	yvalue = (phi_face.yvalue(:,2:Ny+1)-phi_face.yvalue(:,1:Ny))./DY;
	zvalue=[];
	cellGrad=FaceVariable(phi.domain, xvalue, yvalue, zvalue);
elseif (d==2.8)
	Nr = phi.domain.dims(1);
	Ntheta = phi.domain.dims(2);
	DR = repmat(phi.domain.cellsize.x(2:end-1), 1, Ntheta);
	DTHETA = repmat(phi.domain.cellsize.y(2:end-1)', Nr, 1);
	rp = repmat(phi.domain.cellcenters.x, 1, Ntheta);
	xvalue = (phi_face.xvalue(2:Nr+1,:)-phi_face.xvalue(1:Nr,:))./DR;
	yvalue = (phi_face.yvalue(:,2:Ntheta+1)-phi_face.yvalue(:,1:Ntheta))./(DTHETA.*rp);
	zvalue=[];
	cellGrad=FaceVariable(phi.domain, xvalue, yvalue, zvalue);
elseif d == 3
	Nx = phi.domain.dims(1);
	Ny = phi.domain.dims(2);
	Nz = phi.domain.dims(3);
	DX = repmat(phi.domain.cellsize.x(2:end-1), 1, Ny, Nz);
	DY = repmat(phi.domain.cellsize.y(2:end-1)', Nx, 1, Nz);
	DZ = zeros(1,1,Nz);
	DZ(1,1,:) = phi.domain.cellsize.z(2:end-1);
	DZ=repmat(DZ, Nx, Ny, 1);
	xvalue = (phi_face.xvalue(2:Nx+1,:,:)-phi_face.xvalue(1:Nx,:,:))./DX;
	yvalue = (phi_face.yvalue(:,2:Ny+1,:)-phi_face.yvalue(:,1:Ny,:))./DY;
	zvalue = (phi_face.zvalue(:,:,2:Nz+1)-phi_face.zvalue(:,:,1:Nz))./DZ;
	cellGrad=FaceVariable(phi.domain, xvalue, yvalue, zvalue);
elseif d == 3.2
	Nr = phi.domain.dims(1);
	Ntheta = phi.domain.dims(2);
	Nz = phi.domain.dims(3);
	DR = repmat(phi.domain.cellsize.x(2:end-1), 1, Ntheta, Nz);
	DTHETA = repmat(phi.domain.cellsize.y(2:end-1)', Nr, 1, Nz);
	DZ = zeros(1,1,Nz);
	DZ(1,1,:) = phi.domain.cellsize.z(2:end-1);
	DZ = repmat(DZ, Nr, Ntheta, 1);
	rp = repmat(phi.domain.cellcenters.x, 1, Ntheta, Nz);
	xvalue = (phi_face.xvalue(2:Nr+1,:,:)-phi_face.xvalue(1:Nr,:,:))./DR;
	yvalue = (phi_face.yvalue(:,2:Ntheta+1,:)-phi_face.yvalue(:,1:Ntheta,:))./(DTHETA.*rp);
	zvalue = (phi_face.zvalue(:,:,2:Nz+1)-phi_face.zvalue(:,:,1:Nz))./DZ;
	cellGrad=FaceVariable(phi.domain, xvalue, yvalue, zvalue);
end
