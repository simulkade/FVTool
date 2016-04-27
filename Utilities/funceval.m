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

%{
Copyright (c) 2012, 2013, Ali Akbar Eftekhari
All rights reserved.

Redistribution and use in source and binary forms, with or
without modification, are permitted provided that the following
conditions are met:

    *   Redistributions of source code must retain the above copyright notice,
        this list of conditions and the following disclaimer.
    *   Redistributions in binary form must reproduce the above
        copyright notice, this list of conditions and the following
        disclaimer in the documentation and/or other materials provided
        with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR
CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA,
OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
%}

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
end
