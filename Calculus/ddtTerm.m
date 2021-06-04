function ddt = ddtTerm(phi, phi_old, dt, BC)
% function ddt = ddtTerm(dt, phi)
% it returns the derivative of phi with respect to time in the vector ddt
% phi is a structure with elements phi.value and phi.Old
% mesh structure is used to map the (phi.value-phi.Old) matrix in the
% ddt vector
%
% SYNOPSIS:
%		ddt = ddtTerm(phi, phi_old, dt, BC)
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


d = phi.domain.dimension;
if (d ==1) || (d==1.5) || (d==1.8)
	ddt = createCellVariable(phi.domain, (phi.value(2:end-1)-phi_old.value(2:end-1))/dt, BC);
elseif (d == 2) || (d == 2.5) || (d==2.8)
	ddt = createCellVariable(phi.domain, (phi.value(2:end-1, 2:end-1)-phi_old.value(2:end-1, 2:end-1))/dt, BC);
elseif (d == 3) || (d==3.2)
    ddt = createCellVariable(phi.domain, (phi.value(2:end-1, 2:end-1, 2:end-1)-phi_old.value(2:end-1, 2:end-1, 2:end-1))/dt, BC);
end
