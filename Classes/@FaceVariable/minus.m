function r = minus(p,q)
%MINUS: this function subtracts the x, y, and z values of the structures that I use in
% the FVtool.
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

% First implementation:
% m = size(p.xvalue);
% d = ndims(p.xvalue);
% if (min(m) == 1) % 1D: only x value
%     r.xvalue = p.xvalue-q.xvalue;
% elseif (d == 2) % 2D: x and y values
%     r.xvalue = p.xvalue-q.xvalue;
%     r.yvalue = p.yvalue-q.yvalue;
% else % 3D:x, y, and z values
%     r.xvalue = p.xvalue-q.xvalue;
%     r.yvalue = p.yvalue-q.yvalue;
%     r.zvalue = p.zvalue-q.zvalue;
% end
r = p+(-q);
