function [RHSdiv, RHSdivx, RHSdivy] = ...
            divergenceTermCylindrical2D(F)
% This function calculates the divergence of a field using its face
% average value and the vector u, which is a face vector
%
% SYNOPSIS:
%   [RHSdiv, RHSdivx, RHSdivy] = divergenceTermCylindrical2D(F)
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
Nr = F.domain.dims(1);
Nz = F.domain.dims(2);
G=reshape(1:(Nr+2)*(Nz+2), Nr+2, Nz+2);
dr = repmat(F.domain.cellsize.x(2:end-1), 1, Nz);
dz = repmat(F.domain.cellsize.y(2:end-1)', Nr, 1);
rp = repmat(F.domain.cellcenters.x, 1, Nz);
rf = repmat(F.domain.facecenters.x, 1, Nz);

% define the vector of cell index
row_index = reshape(G(2:Nr+1,2:Nz+1),Nr*Nz,1); % main diagonal (only internal cells)

% calculate the flux vector
% note: size(Fx) = [1:m+1, 1:n] and size(Fy) = [1:m, 1:n+1]
Fx = F.xvalue;
Fy = F.yvalue;

% reassign the east, west, north, and south flux vectors for the
% code readability
Fe = Fx(2:Nr+1,:);		Fw = Fx(1:Nr,:);
Fn = Fy(:,2:Nz+1);       Fs = Fy(:,1:Nz);
re = rf(2:Nr+1,:);         rw = rf(1:Nr,:);

% compute the divergence
div_x = (re.*Fe-rw.*Fw)./(rp.*dr);
div_y = (Fn-Fs)./dz;

% define the RHS Vector
RHSdiv = zeros((Nr+2)*(Nz+2),1);
RHSdivx = zeros((Nr+2)*(Nz+2),1);
RHSdivy = zeros((Nr+2)*(Nz+2),1);

% assign the values of the RHS vector
RHSdiv(row_index) = reshape(div_x+div_y,Nr*Nz,1);
RHSdivx(row_index) = reshape(div_x,Nr*Nz,1);
RHSdivy(row_index) = reshape(div_y,Nr*Nz,1);
