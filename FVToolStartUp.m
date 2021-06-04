function FVToolStartUp()
%
% SYNOPSIS:
%   FVToolStartUp()
%
% PARAMETERS:
%   No perameter
%
% RETURNS:
%   None
%
% EXAMPLE:
%   n.a.
%
% SEE ALSO:
%   PVTinitialize, FVTdemo

%{
Copyright (c) 2012-2021 Ali Akbar Eftekhari
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
try
p = mfilename('fullpath');
file_name = mfilename;
current_path = p(1:end-1-length(file_name));
addpath([current_path '/Boundary']);
addpath([current_path '/Calculus']);
addpath([current_path '/Classes']);
addpath([current_path '/Discretization']);
addpath([current_path '/MeshGeneration']);
addpath([current_path '/Solvers']);
addpath([current_path '/Utilities']);
addpath([current_path '/Visualization']);
addpath([current_path '/Examples']);
addpath([current_path '/Physics']);
addpath([current_path '/Tests']);

try
    addpath([current_path '/PhysicalProperties']);
    addpath([current_path '/FieldGeology']);
catch
    disp(['Some of the physical functions are not available in this copy.' ...
    ' It does not affect the functionality of the FVMtool']);
end

% Check for other solvers
    % check for AGMG availability
    cd(current_path);
    cd('Solvers');
    if exist('AGMG_3.2', 'dir') == 7
        addpath([pwd '/AGMG_3.2']);
        disp('AGMG 3.2 linear solver is available.');
    elseif exist('AGMG_3.0', 'dir') == 7
        addpath([pwd '/AGMG_3.0']);
        disp('AGMG 3.0 linear solver is available.');
    else
        disp('AGMG 3.x linear solver is NOT available (Not necessary).');
    end
    % check for Factorize availability
    if exist('Factorize', 'dir') == 7
        addpath([pwd '/Factorize']);
        disp('Factorize is available.');
    end
    cd(current_path);
% end of check for other solvers

% Check for the PVT package
if exist('PVTtoolbox', 'dir') == 7
    addpath([current_path '/PVTtoolbox']);
    if exist('PVTinitialize', 'file') == 2
        cd('PVTtoolbox');
        PVTinitialize();
        cd(current_path);
        disp('PVTtoolbox has started successfully.');
    else
        disp('PVTtoolbox is found but cannot be initialized.');
    end
else
    disp('PVTtoolbox is NOT available (Not necessary).');
end
disp('FiniteVolumeToolbox has started successfully.');
catch err
    error('An error occured while tryng to start the FiniteVolumeToolbox.')
end
