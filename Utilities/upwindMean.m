function phiFaceAverage = upwindMean(phi, u)
% This function gets the value of the field variable phi defined
% over the MeshStructure and calculates the upwind average on
% the cell faces, for a uniform mesh based on the direction of the velocity
% vector u
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

d = phi.domain.dimension;
if (d ==1) || (d==1.5) || (d==1.8)
	phiFaceAverage = upwindMean1D(phi, u);
elseif (d == 2) || (d == 2.5) || (d==2.8)
	phiFaceAverage = upwindMean2D(phi, u);
elseif (d == 3) || (d==3.2)
    phiFaceAverage = upwindMean3D(phi, u);
end
