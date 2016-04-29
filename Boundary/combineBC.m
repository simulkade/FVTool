function [Mout, RHSout] = combineBC(MeshStructure, BC, Meq, RHSeq)
%COMBINEBC This function combines the boundary condition equations with the
%main physical model equations, and delivers the matrix of coefficient and
%RHS to be solved for the internal cells.
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

d = MeshStructure.dimension;

if (d ==1) || (d==1.5)
	[Mout, RHSout] = combineBC1D(MeshStructure, BC, Meq, RHSeq);
elseif (d == 2) || (d == 2.5)
	[Mout, RHSout] = combineBC2D(MeshStructure, BC, Meq, RHSeq);
elseif d == 3
    [Mout, RHSout] = combineBC3D(MeshStructure, BC, Meq, RHSeq);
end
