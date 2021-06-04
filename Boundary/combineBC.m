function [Mout, RHSout] = combineBC(BC, Meq, RHSeq)
%COMBINEBC This function combines the boundary condition equations with the
%main physical model equations, and delivers the matrix of coefficient and
%RHS to be solved for the internal cells. It is useful if one needs to use 
%and ODE solver for the accumulation term, i.e.
% d phi/ dt = M phi
%
% SYNOPSIS:
%    [Mout, RHSout] = combineBC(BC, Meq, RHSeq)
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


d = BC.domain.dimension;

if (d ==1) || (d==1.5)
	[Mout, RHSout] = combineBC1D(BC, Meq, RHSeq);
elseif (d == 2) || (d == 2.5)
	[Mout, RHSout] = combineBC2D(BC, Meq, RHSeq);
elseif d == 3
    [Mout, RHSout] = combineBC3D(BC, Meq, RHSeq);
end
