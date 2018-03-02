function RHS = convectionTvdRHSSpherical1D(u, phi, FL)
% This function uses the upwind scheme to discretize a 1D
% convection term in the form \grad (u \phi) where u is a face vactor
% It also returns the x and y parts of the matrix of coefficient.
%
% SYNOPSIS:
%   [M, RHS] = convectionTvdTermCylindrical1D(u, phi, FL)
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

% extract data from the mesh structure
Nr = u.domain.dims(1);
DXp = u.domain.cellsize.x(2:end-1);
dx = 0.5*(u.domain.cellsize.x(1:end-1)+u.domain.cellsize.x(2:end));
r = u.domain.cellcenters.x;
rf = u.domain.facecenters.x;
RHS = zeros(Nr+2, 1);
psi_p = zeros(Nr+1,1);
psi_m = zeros(Nr+1,1);

% extract the velocity data
% note: size(ux) = [1:m+1, 1:n] and size(uy) = [1:m, 1:n+1]
ux = u.xvalue;

% calculate the upstream to downstream gradient ratios for u>0 (+ ratio)
% P is 3:Nr+2
% W is 2:Nr+1
% WW is 1:Nr
dphi_p = (phi.value(2:Nr+2)-phi.value(1:Nr+1))./dx;
rp = dphi_p(1:end-1)./fsign(dphi_p(2:end));
psi_p(2:Nr+1) = 0.5*FL(rp).*(phi.value(3:Nr+2)-phi.value(2:Nr+1));
psi_p(1) = 0.0; % left boundary will be handled explicitly

% calculate the upstream to downstream gradient ratios for u<0 (- ratio)
% P is 3:Nr+2
% W is 2:Nr+1
% WW is 1:Nr
rm = dphi_p(2:end)./fsign(dphi_p(1:end-1));
psi_m(1:Nr) = 0.5*FL(rm).*(phi.value(1:Nr)-phi.value(2:Nr+1));
psi_m(Nr+1) = 0.0; % right boundary will be handled explicitly

% reassign the east, west, north, and south velocity vectors for the
% code readability
ue = ux(2:Nr+1);		uw = ux(1:Nr);
re = rf(2:Nr+1);     rw = rf(1:Nr);

% find the velocity direction for the upwind scheme
ue_min = min(ue,0);	ue_max = max(ue,0);
uw_min = min(uw,0);	uw_max = max(uw,0);

% calculate the TVD correction term
RHS(2:Nr+1) = -(1./(1/3*(rf(2:Nr+1).^3-rf(1:Nr).^3))).*(re.^2.*(ue_max.*psi_p(2:Nr+1)+ue_min.*psi_m(2:Nr+1))- ...
              rw.^2.*(uw_max.*psi_p(1:Nr)+uw_min.*psi_m(1:Nr)));

end

function phi_out = fsign(phi_in)
% This function checks the value of phi_in and assigns an eps value to the
% elements that are less than or equal to zero, while keeping the signs of
% the nonzero elements
    phi_out = (abs(phi_in)>=eps).*phi_in+eps*(phi_in==0)+eps*(abs(phi_in)<eps).*sign(phi_in);
end
