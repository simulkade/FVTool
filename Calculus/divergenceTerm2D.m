function [RHSdiv, RHSdivx, RHSdivy] = divergenceTerm2D(F)
% This function calculates the divergence of a field using its face
% average flux vector facevariable, which is a face vector
%
% SYNOPSIS:
%   [RHSdiv, RHSdivx, RHSdivy] = divergenceTerm2D(F)
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
Nx = F.domain.dims(1);
Ny = F.domain.dims(2);
G=reshape(1:(Nx+2)*(Ny+2), Nx+2, Ny+2);
DX = repmat(F.domain.cellsize.x(2:end-1), 1, Ny);
DY = repmat(F.domain.cellsize.y(2:end-1)', Nx, 1);

% define the vector of cell index
row_index = reshape(G(2:Nx+1,2:Ny+1),Nx*Ny,1); % main diagonal (only internal cells)


% reassign the east, west, north, and south flux vectors for the
% code readability
Fe = F.xvalue(2:Nx+1,:);		Fw = F.xvalue(1:Nx,:);
Fn = F.yvalue(:,2:Ny+1);       Fs = F.yvalue(:,1:Ny);

% compute the divergence
div_x = (Fe - Fw)./DX;
div_y = (Fn - Fs)./DY;

% define the RHS Vector
RHSdiv = zeros((Nx+2)*(Ny+2),1);
RHSdivx = zeros((Nx+2)*(Ny+2),1);
RHSdivy = zeros((Nx+2)*(Ny+2),1);

% assign the values of the RHS vector
RHSdiv(row_index) = reshape(div_x+div_y,Nx*Ny,1);
RHSdivx(row_index) = reshape(div_x,Nx*Ny,1);
RHSdivy(row_index) = reshape(div_y,Nx*Ny,1);
