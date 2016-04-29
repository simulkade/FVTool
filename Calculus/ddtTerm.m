function ddt = ddtTerm(MeshStructure, dt, phi)
% function ddt = ddtTerm(MeshStructure, dt, phi)
% it returns the derivative of phi with respect to time in the vector ddt
% phi is a structure with elements phi.value and phi.Old
% mesh structure is used to map the (phi.value-phi.Old) matrix in the
% ddt vector
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
	ddt = ddtTerm1D(MeshStructure, dt, phi);
elseif (d == 2) || (d == 2.5) || (d==2.8)
	ddt = ddtTerm2D(MeshStructure, dt, phi);
elseif (d == 3) || (d==3.2)
    ddt = ddtTerm3D(MeshStructure, dt, phi);
end
