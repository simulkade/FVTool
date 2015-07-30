function [M, RHS] = transientTerm(phi_old, dt, varargin)
% function [M, RHS] = transientTerm(MeshStructure, h, dt, phi)
% Matrix of coefficients and the RHS vector for a transient term 
% alfa \partial_t \phi
% 
% 
% SYNOPSIS:
%   [M, RHS] = transientTerm(phi_old, dt)
%   [M, RHS] = transientTerm(phi_old, dt, alfa)
% 
% PARAMETERS:
%   phi_old: Cell Variable
%   dt:      time step
%   alfa:    Cell Variable
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


d = phi_old.domain.dimension;
if (d ==1) || (d==1.5)
	% extract data from the mesh structure
    Nx = phi_old.domain.dims(1);
    G = 1:Nx+2;
    if nargin==2
        alfa=ones(Nx,1);
    elseif isa(varargin{1}, 'CellVariable')
        alfa=varargin{1}.value(2:Nx+1);
    else
        alfa=varargin{1}.*ones(Nx,1);
    end
    % rearrange the matrix of k and build the sparse matrix for internal cells
    row_index = reshape(G(2:Nx+1),Nx,1); % main diagonal (only internal cells)
    AP_diag = reshape(alfa/dt,Nx,1);
    M = sparse(row_index, row_index, AP_diag, Nx+2, Nx+2);

    % define the RHS Vector
    RHS = zeros(Nx+2,1);

    % assign the values of the RHS vector
    RHS(row_index) = reshape(alfa.*phi_old.value(2:Nx+1)/dt,Nx,1);
elseif (d == 2) || (d == 2.5) || (d == 2.8)
	Nxy = phi_old.domain.dims;
    Nx = Nxy(1); Ny = Nxy(2);
    G=reshape(1:(Nx+2)*(Ny+2), Nx+2, Ny+2);
    if nargin==2
        alfa=ones(Nx,Ny);
    elseif isa(varargin{1}, 'CellVariable')
        alfa=varargin{1}.value(2:Nx+1, 2:Ny+1);
    else
        alfa=varargin{1}.*ones(Nx,Ny);
    end

    % rearrange the matrix of k and build the sparse matrix for internal cells
    row_index = reshape(G(2:Nx+1,2:Ny+1),Nx*Ny,1); % main diagonal (only internal cells)
    AP_diag = reshape(alfa/dt,Nx*Ny,1);
    M = sparse(row_index, row_index, AP_diag, (Nx+2)*(Ny+2), (Nx+2)*(Ny+2));

    % define the RHS Vector
    RHS = zeros((Nx+2)*(Ny+2),1);

    % assign the values of the RHS vector
    RHS(row_index) = reshape(alfa.*phi_old.value(2:Nx+1,2:Ny+1)/dt,Nx*Ny,1);
elseif (d == 3) || (d == 3.2)
    Nxyz = phi_old.domain.dims;
    Nx = Nxyz(1); Ny = Nxyz(2); Nz = Nxyz(3);
    G=reshape(1:(Nx+2)*(Ny+2)*(Nz+2), Nx+2, Ny+2, Nz+2);
    if nargin==2
        alfa=ones(Nx,Ny,Nz);
    elseif isa(varargin{1}, 'CellVariable')
        alfa=varargin{1}.value(2:Nx+1, 2:Ny+1, 2:Nz+1);
    else
        alfa=varargin{1}.*ones(Nx,Ny,Nz);
    end
    % rearrange the matrix of k and build the sparse matrix for internal cells
    row_index = reshape(G(2:Nx+1,2:Ny+1,2:Nz+1),Nx*Ny*Nz,1); % main diagonal (only internal cells)
    AP_diag = reshape(alfa/dt,Nx*Ny*Nz,1);
    M = sparse(row_index, row_index, AP_diag, ...
        (Nx+2)*(Ny+2)*(Nz+2), (Nx+2)*(Ny+2)*(Nz+2));

    % define the RHS Vector
    RHS = zeros((Nx+2)*(Ny+2)*(Nz+2),1);

    % assign the values of the RHS vector
    RHS(row_index) = reshape(alfa.*phi_old.value(2:Nx+1,2:Ny+1,2:Nz+1)/dt,Nx*Ny*Nz,1);
end
