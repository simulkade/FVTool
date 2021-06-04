function visualizeCellVectors(phi_cell)
%VISUALIZECELLS plots the values of cell variable phi.value
%
% SYNOPSIS:
%  visualizeCellVectors(phi_cell)
%
% PARAMETERS:
%  phi_cell: a CellVector variable
%
% RETURNS:
%  None
%
% EXAMPLE:
%
% SEE ALSO:
%

% Written by Ali A. Eftekhari
% See the license file

d = phi_cell.domain.dimension;
switch d
    case 1
        warning('No vector visualization for a 1D domain.');
    case 1.5
        warning('No vector visualization for a 1D domain.');
    case 2
        visualizeCellVectors2D(phi_cell);
    case 2.5
        visualizeCellVectors2D(phi_cell);
    case 2.8
        visualizeCellVectorsRadial2D(phi_cell);
    case 3
        visualizeCellVectors3D(phi_cell);
    case 3.2
        visualizeCellVectorsCylindrical3D(phi_cell);
end

end
