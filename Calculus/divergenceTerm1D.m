function RHSdiv = divergenceTerm1D(F)
% This function calculates the divergence of a field using its face
% average value
%
% SYNOPSIS:
%   RHSdiv = divergenceTerm1D(F)
%
% PARAMETERS:
%   F: Face Variable
%
% RETURNS:
%
%
% EXAMPLE:
%
% SEE ALSO:
%


% extract data from the mesh structure
Nx = F.domain.dims(1);
G = 1:Nx+2;
DX = F.domain.cellsize.x(2:end-1);

% define the vector of cell index
row_index = reshape(G(2:Nx+1),Nx,1); % main diagonal (only internal cells)

% compute the divergence
div_x = (F.xvalue(2:Nx+1)-F.xvalue(1:Nx))./DX;

% define the RHS Vector
RHSdiv = zeros(Nx+2,1);

% assign the values of the RHS vector
RHSdiv(row_index) = reshape(div_x,Nx,1);
