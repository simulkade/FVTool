function [RHSdiv, RHSdivx, RHSdivy] = ...
            divergenceTermRadial2D(F)
% This function calculates the divergence of a field using its face
% average value and the vector u, which is a face vector
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

% extract data from the mesh structure
Nr = F.domain.dims(1);
Ntetta = F.domain.dims(2);
G=reshape(1:(Nr+2)*(Ntetta+2), Nr+2, Ntetta+2);
dr = repmat(F.domain.cellsize.x(2:end-1), 1, Ntetta);
dtetta = repmat(F.domain.cellsize.y(2:end-1)', Nr, 1);
rp = repmat(F.domain.cellcenters.x, 1, Ntetta);
rf = repmat(F.domain.facecenters.x, 1, Ntetta);
re = rf(2:Nr+1,:);         rw = rf(1:Nr,:);

% define the vector of cell index
row_index = reshape(G(2:Nr+1,2:Ntetta+1),Nr*Ntetta,1); % main diagonal (only internal cells)

% reassign the east, west, north, and south flux vectors for the
% code readability
Fe = F.xvalue(2:Nr+1,:);		Fw = F.xvalue(1:Nr,:);
Fn = F.yvalue(:,2:Ntetta+1);       Fs = F.yvalue(:,1:Ntetta);

% compute the divergence
div_x = (re.*Fe-rw.*Fw)./(rp.*dr);
div_y = (Fn-Fs)./(dtetta.*rp);

% define the RHS Vector
RHSdiv = zeros((Nr+2)*(Ntetta+2),1);
RHSdivx = zeros((Nr+2)*(Ntetta+2),1);
RHSdivy = zeros((Nr+2)*(Ntetta+2),1);

% assign the values of the RHS vector
RHSdiv(row_index) = reshape(div_x+div_y,Nr*Ntetta,1);
RHSdivx(row_index) = reshape(div_x,Nr*Ntetta,1);
RHSdivy(row_index) = reshape(div_y,Nr*Ntetta,1);
