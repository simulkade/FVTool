function r = funceval(f, varargin)
%FUNCEVAL calculates the value of the function f at the face variable phi
% f must be able to accept matrices
% f: function handle
% phi: input to the function f
% r: output of the function f
% f can handle 3 input phi variables, but they should all be of the same
% class
%
% SYNOPSIS:
%   r = funceval(f, varargin)
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

% Written by Ali A. Eftekhari
% See the license file

% analyze phi to find the dimension of the matrix

switch nargin
    case 2
        r= FaceVariable(varargin{1}.domain, ...
        f(varargin{1}.xvalue), ...
        f(varargin{1}.yvalue), ...
        f(varargin{1}.zvalue));
    case 3
        r= FaceVariable(varargin{1}.domain, ...
        f(varargin{1}.xvalue, varargin{2}.xvalue), ...
        f(varargin{1}.yvalue, varargin{2}.yvalue), ...
        f(varargin{1}.zvalue, varargin{2}.zvalue));
    case 4
        r= FaceVariable(varargin{1}.domain, ...
        f(varargin{1}.xvalue, varargin{2}.xvalue, varargin{3}.xvalue), ...
        f(varargin{1}.yvalue, varargin{2}.yvalue, varargin{3}.yvalue), ...
        f(varargin{1}.zvalue, varargin{2}.zvalue, varargin{3}.zvalue));
    case 5
        r= FaceVariable(varargin{1}.domain, ...
        f(varargin{1}.xvalue, varargin{2}.xvalue, varargin{3}.xvalue, varargin{4}.xvalue), ...
        f(varargin{1}.yvalue, varargin{2}.yvalue, varargin{3}.yvalue, varargin{4}.yvalue), ...
        f(varargin{1}.zvalue, varargin{2}.zvalue, varargin{3}.zvalue, varargin{4}.zvalue));
    case 6
        r= FaceVariable(varargin{1}.domain, ...
        f(varargin{1}.xvalue, varargin{2}.xvalue, varargin{3}.xvalue, varargin{4}.xvalue, varargin{5}.xvalue),...
        f(varargin{1}.yvalue, varargin{2}.yvalue, varargin{3}.yvalue, varargin{4}.yvalue, varargin{5}.yvalue),...
        f(varargin{1}.zvalue, varargin{2}.zvalue, varargin{3}.zvalue, varargin{4}.zvalue, varargin{5}.zvalue));
end
