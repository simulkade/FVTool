function r = plus(p,q)
% this function adds the x, y, and z values of the structures that I use in
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

if (isa(p, 'FaceVariable')&&isa(q, 'FaceVariable'))
    r=p;
    r.xvalue = p.xvalue+q.xvalue;
    r.yvalue = p.yvalue+q.yvalue;
    r.zvalue = p.zvalue+q.zvalue;
elseif isa(p, 'FaceVariable')
    r=p;
    r.xvalue = p.xvalue+q;
    r.yvalue = p.yvalue+q;
    r.zvalue = p.zvalue+q;
else
    r=q;
    r.xvalue = p+q.xvalue;
    r.yvalue = p+q.yvalue;
    r.zvalue = p+q.zvalue;
end
