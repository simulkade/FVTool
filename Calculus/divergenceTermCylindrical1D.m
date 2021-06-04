function RHSdiv = divergenceTermCylindrical1D(F)
% This function calculates the divergence of a field using its face
% average value and the vector u, which is a face vector
%
% SYNOPSIS:
%       RHSdiv = divergenceTermCylindrical1D(F)
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

% extract data from the mesh structure
Nr = F.domain.dims(1);
G = 1:Nr+2;
DX = F.domain.cellsize.x(2:end-1);
rp = F.domain.cellcenters.x;
rf = F.domain.facecenters.x;

% define the vector of cell index
row_index = reshape(G(2:Nr+1),Nr,1); % main diagonal (only internal cells)

% calculate the flux vector
% note: size(Fx) = [1:m+1]
Fx = F.xvalue;

% reassign the east, west, north, and south flux vectors for the
% code readability
Fe = Fx(2:Nr+1);		Fw = Fx(1:Nr);
re = rf(2:Nr+1);     rw = rf(1:Nr);

% compute the divergence
div_x = (re.*Fe - rw.*Fw)./(rp.*DX);

% define the RHS Vector
RHSdiv = zeros(Nr+2,1);

% assign the values of the RHS vector
RHSdiv(row_index) = reshape(div_x,Nr,1);
