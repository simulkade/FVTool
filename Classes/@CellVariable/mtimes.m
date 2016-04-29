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

if (isa(p, 'CellVariable')&&isa(q, 'CellVariable'))
    error('FVMtool: Wrong use of mtimes for a cell variable. Try using .* instead.');
elseif isa(p, 'CellVariable')
    r=p;
    r.value = p.value*q;
else
    r=q;
    r.value = p*q.value;
end
