function r = mtimes(p,q)
%TIMES this function multiplies the x, y, and z values of the structures that I use in
% the FVtool.
% reassign the scalar to new variable a and the struct to new variable b
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

if (isa(p, 'CellVector')&&isa(q, 'CellVector'))
    error('FVMtool: Wrong use of mtimes for a cell vector. Try using .* instead.');
elseif isa(p, 'CellVector')
    r=p;
    r.xvalue = p.xvalue*q;
    r.yvalue = p.yvalue*q;
    r.zvalue = p.zvalue*q;
else
    r=q;
    r.xvalue = p*q.xvalue;
    r.yvalue = p*q.yvalue;
    r.zvalue = p*q.zvalue;
end
