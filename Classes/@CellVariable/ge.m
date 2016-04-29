function r = ge(p,q)
% this function compares the greater than or equality condition of x, y, and z values of the structures that I use in
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

if (isa(p, 'CellVariable')&&isa(q, 'CellVariable'))
    r=p;
    r.value = p.value>=q.value;
elseif isa(p, 'CellVariable')
    r=p;
    r.value = p.value>=q;
else
    r=q;
    r.value = p>=q.value;
end
