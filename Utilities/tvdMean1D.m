function phiFaceAverage = tvdMean1D(MeshStructure, phi, u, FL)
% This function gets the value of the field variable phi defined
% over the MeshStructure and calculates the TVD average on 
% the cell faces, based on the direction of the velocity vector for a uniform mesh.
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
Nx = MeshStructure.numberofcells;

% extract the velocity data
% note: size(ux) = [1:m+1, 1:n] and size(uy) = [1:m, 1:n+1]
ux = u.xvalue;
phi_p = zeros(Nx+1,1);
phi_m = zeros(Nx+1,1);

% calculate the upstream to downstream gradient ratios for u>0 (+ ratio)
% P is 3:Nx+2
% W is 2:Nx+1
% WW is 1:Nx
dphi_p = phi(2:Nx+2)-phi(1:Nx+1);
rp = dphi_p(1:end-1)./fsign(dphi_p(2:end));
phi_p(2:Nx+1) = phi(2:Nx+1)+0.5*FL(rp).*(phi(3:Nx+2)-phi(2:Nx+1));
phi_p(1) = (phi(1)+phi(2))/2; % left boundary

% calculate the upstream to downstream gradient ratios for u<0 (- ratio)
% P is 3:Nx+2
% W is 2:Nx+1
% WW is 1:Nx
rm = dphi_p(2:end)./fsign(dphi_p(1:end-1));
phi_m(1:Nx) = phi(2:Nx+1)+0.5*FL(rm).*(phi(1:Nx)-phi(2:Nx+1));
phi_m(Nx+1) = (phi(end)+phi(end-1))/2; % right boundary

% calculate the average value
phiFaceAverage.xvalue = (ux>0).*phi_p+ ...
                        (ux<0).*phi_m+ ...
                        0.5*(ux==0).*(phi(1:Nx+1)+phi(2:Nx+2));
end

function phi_out = fsign(phi_in)
% This function checks the value of phi_in and assigns an eps value to the
% elements that are less than or equal to zero, while keeping the signs of
% the nonzero elements
    phi_out = (abs(phi_in)>=eps).*phi_in+eps*(phi_in==0)+eps*(abs(phi_in)<eps).*sign(phi_in);
end
