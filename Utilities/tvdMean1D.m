function phiFaceAverage = tvdMean1D(phi, u, FL)
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

% Copyright (c) 2012-2016 Ali Akbar Eftekhari
% See the license file

% check the size of the variable and the mesh dimension
Nx = phi.domain.dims(1);
dx = 0.5*(u.domain.cellsize.x(1:end-1)+u.domain.cellsize.x(2:end));

% extract the velocity data
% note: size(ux) = [1:m+1, 1:n] and size(uy) = [1:m, 1:n+1]
ux = u.xvalue;
phi_p = zeros(Nx+1,1);
phi_m = zeros(Nx+1,1);

% calculate the upstream to downstream gradient ratios for u>0 (+ ratio)
% P is 3:Nx+2
% W is 2:Nx+1
% WW is 1:Nx
dphi_p = (phi.value(2:Nx+2)-phi.value(1:Nx+1))./dx;
rp = dphi_p(1:end-1)./fsign(dphi_p(2:end));
phi_p(2:Nx+1) = phi.value(2:Nx+1)+0.5*FL(rp).*(phi.value(3:Nx+2)-phi.value(2:Nx+1));
phi_p(1) = (phi.value(1)+phi.value(2))/2; % left boundary

% calculate the upstream to downstream gradient ratios for u<0 (- ratio)
% P is 3:Nx+2
% W is 2:Nx+1
% WW is 1:Nx
rm = dphi_p(2:end)./fsign(dphi_p(1:end-1));
phi_m(1:Nx) = phi.value(2:Nx+1)+0.5*FL(rm).*(phi.value(1:Nx)-phi.value(2:Nx+1));
phi_m(Nx+1) = (phi.value(end)+phi.value(end-1))/2; % right boundary

% calculate the average value
xvalue = (ux>0).*phi_p+ ...
                        (ux<0).*phi_m+ ...
                        0.5*(ux==0).*(phi.value(1:Nx+1)+phi.value(2:Nx+2));

phiFaceAverage=FaceVariable(phi.domain, xvalue, [],[]);
end

function phi_out = fsign(phi_in)
% This function checks the value of phi_in and assigns an eps value to the
% elements that are less than or equal to zero, while keeping the signs of
% the nonzero elements
    phi_out = (abs(phi_in)>=eps).*phi_in+eps*(phi_in==0)+eps*(abs(phi_in)<eps).*sign(phi_in);
end
