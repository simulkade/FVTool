function phiFaceAverage = tvdMean2D(phi, u, FL)
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

% extract the velocity data
% note: size(ux) = [1:m+1, 1:n] and size(uy) = [1:m, 1:n+1]
ux = u.xvalue;
uy = u.yvalue;

% check the size of the variable and the mesh dimension
Nx = u.domain.dims(1);
Ny = u.domain.dims(2);
dx=repmat(0.5*(u.domain.cellsize.x(1:end-1)+u.domain.cellsize.x(2:end)), 1, Ny);
dy=repmat(0.5*(u.domain.cellsize.y(1:end-1)+u.domain.cellsize.y(2:end))', Nx, 1);

%
phiX_p = zeros(Nx+1,Ny);
phiX_m = zeros(Nx+1,Ny);
phiY_p = zeros(Nx,Ny+1);
phiY_m = zeros(Nx,Ny+1);

% calculate the upstream to downstream gradient ratios for u>0 (+ ratio)
% x direction
dphiX_p = (phi.value(2:Nx+2, 2:Ny+1)-phi.value(1:Nx+1, 2:Ny+1))./dx;
rX_p = dphiX_p(1:end-1,:)./fsign(dphiX_p(2:end,:));
phiX_p(2:Nx+1,:) = phi.value(2:Nx+1, 2:Ny+1)+0.5*FL(rX_p).* ...
    (phi.value(3:Nx+2,2:Ny+1)-phi.value(2:Nx+1, 2:Ny+1));
phiX_p(1, :) = (phi.value(1, 2:Ny+1)+phi.value(2, 2:Ny+1))/2; % left boundary
% y direction
dphiY_p = (phi.value(2:Nx+1, 2:Ny+2)-phi.value(2:Nx+1, 1:Ny+1))./dy;
rY_p = dphiY_p(:,1:end-1)./fsign(dphiY_p(:,2:end));
phiY_p(:,2:Ny+1) = phi.value(2:Nx+1, 2:Ny+1)+0.5*FL(rY_p).* ...
    (phi.value(2:Nx+1,3:Ny+2)-phi.value(2:Nx+1, 2:Ny+1));
phiY_p(:,1) = (phi.value(2:Nx+1,1)+phi.value(2:Nx+1,2))/2; % Bottom boundary

% calculate the upstream to downstream gradient ratios for u<0 (- ratio)
% x direction
rX_m = dphiX_p(2:end,:)./fsign(dphiX_p(1:end-1,:));
phiX_m(1:Nx,:) = phi.value(2:Nx+1, 2:Ny+1)+0.5*FL(rX_m).* ...
    (phi.value(1:Nx, 2:Ny+1)-phi.value(2:Nx+1, 2:Ny+1));
phiX_m(Nx+1,:) = (phi.value(end, 2:Ny+1)+phi.value(end-1, 2:Ny+1))/2; % right boundary
% y direction
rY_m = dphiY_p(:,2:end)./fsign(dphiY_p(:,1:end-1));
phiY_m(:,1:Ny) = phi.value(2:Nx+1, 2:Ny+1)+0.5*FL(rY_m).* ...
    (phi.value(2:Nx+1, 1:Ny)-phi.value(2:Nx+1, 2:Ny+1));
phiY_m(:, Ny+1) = (phi.value(2:Nx+1, end)+phi.value(2:Nx+1, end-1))/2; % top boundary


% calculate the average value
xvalue = (ux>0).*phiX_p+ ...
                        (ux<0).*phiX_m+ ...
                        0.5*(ux==0).*(phi.value(1:Nx+1,2:Ny+1)+phi.value(2:Nx+2,2:Ny+1));
yvalue = (uy>0).*phiY_p+ ...
                        (uy<0).*phiY_m+ ...
                        0.5*(uy==0).*(phi.value(2:Nx+1,1:Ny+1)+phi.value(2:Nx+1,2:Ny+2));
phiFaceAverage=FaceVariable(phi.domain, xvalue, yvalue, []);
end

function phi_out = fsign(phi_in)
% This function checks the value of phi_in and assigns an eps value to the
% elements that are less than or equal to zero, while keeping the signs of
% the nonzero elements
    phi_out = (abs(phi_in)>=eps).*phi_in+eps*(phi_in==0)+eps*(abs(phi_in)<eps).*sign(phi_in);
end
