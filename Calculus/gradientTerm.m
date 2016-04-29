function faceGrad = gradientTerm(phi)
% this function calculates the gradient of a variable in x direction
% it checks for the availability of the ghost variables and use them, otherwise
% estimate them, assuming a zero gradient on the boundaries
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

d = phi.domain.dimension;
if (d ==1) || (d==1.5)
	faceGrad = gradientTerm1D(phi);
elseif (d == 2) || (d == 2.5)
	faceGrad = gradientTerm2D(phi);
elseif (d==2.8)
    faceGrad = gradientTermRadial2D(phi);
elseif d == 3
    faceGrad = gradientTerm3D(phi);
elseif d == 3.2
    faceGrad = gradientTermCylindrical3D(phi);
end
