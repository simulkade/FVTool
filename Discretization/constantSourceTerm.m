function RHS = constantSourceTerm(phi)
% RHS vector for an explicit source term
% k is a cell variable
%
% SYNOPSIS:
%   RHS = constantSourceTerm(k)
%
% PARAMETERS:
%   phi: Cell Variable
%
% RETURNS:
%   RHS: vector
%
% EXAMPLE:
%
% SEE ALSO:
%

d = phi.domain.dimension;
if (d ==1) || (d==1.5) || (d==1.8)
	RHS = constantSourceTerm1D(phi);
elseif (d == 2) || (d == 2.5) || (d==2.8)
	RHS = constantSourceTerm2D(phi);
elseif (d==3) || (d==3.2)
    RHS = constantSourceTerm3D(phi);
end
