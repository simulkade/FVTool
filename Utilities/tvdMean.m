function phiFaceAverage = tvdMean(phi, u, FL)
% This function gets the value of the field variable phi defined
% over the MeshStructure and calculates the tvd average on
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
if (d ==1) || (d==1.5)
	phiFaceAverage = tvdMean1D(phi, u, FL);
elseif (d == 2) || (d == 2.5) || (d==2.8)
	phiFaceAverage = tvdMean2D(phi, u, FL);
elseif (d == 3) || (d==3.2)
    phiFaceAverage = tvdMean3D(phi, u, FL);
end
