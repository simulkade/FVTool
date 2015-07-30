% FVTool test script
% This script is supposed to call all the functions of the FVTool package
clc; clear;
%% Part I: creating an array of different mesh types:
% domain size
Lx= 1.0;
Ly= 2*pi;
Lz= 2.0;
Nx=5;
Ny=7;
Nz=9;
X=[0.0 0.1 0.3 0.5 0.55 1.0];
Y= [0.0 0.1 1.0 1.5 2.9 3.0 pi 2*pi];
Z= [0.0 0.01 0.1 0.5 0.7 0.95 1.0 1.25 1.39 2.0];
N_mesh=7;
% create uniform mesh
mesh_uniform= cell(N_mesh,1);
mesh_uniform{1}=createMesh1D(Nx, Lx);
mesh_uniform{2}=createMesh2D(Nx, Ny, Lx, Ly);
mesh_uniform{3}=createMesh3D(Nx, Ny, Nz, Lx, Ly, Lz);
mesh_uniform{4}=createMeshCylindrical1D(Nx, Lx);
mesh_uniform{5}=createMeshCylindrical2D(Nx, Ny, Lx, Ly);
mesh_uniform{6}=createMeshCylindrical3D(Nx, Ny, Nz, Lx, Ly, Lz);
mesh_uniform{7}=createMeshRadial2D(Nx, Ny, Lx, Ly);
disp('Uniform mesh created successfully!');
% create nonuniform mesh
mesh_nonuniform= cell(N_mesh,1);
mesh_nonuniform{1}=createMesh1D(X);
mesh_nonuniform{2}=createMesh2D(X, Y);
mesh_nonuniform{3}=createMesh3D(X, Y, Z);
mesh_nonuniform{4}=createMeshCylindrical1D(X);
mesh_nonuniform{5}=createMeshCylindrical2D(X, Y);
mesh_nonuniform{6}=createMeshCylindrical3D(X, Y, Z);
mesh_nonuniform{7}=createMeshRadial2D(X, Y);
disp('Non-uniform mesh created successfully!');
%% Part II: create cell and face variables
c_val= 0.5;
% uniform
c_u=cell(N_mesh, 1);
for i=1:N_mesh
    c_u{i}= createCellVariable(mesh_uniform{i}, c_val);
end
disp('Cells of fixed values over uniform mesh created successfully!');
for i=1:N_mesh
    c_rand = rand(mesh_uniform{i}.dims);
    c_u{i}= createCellVariable(mesh_uniform{i}, c_rand);
end
disp('Cells of random values over uniform mesh created successfully!');
% nonuniform
c_n=cell(N_mesh, 1);
for i=1:N_mesh
    c_n{i}= createCellVariable(mesh_nonuniform{i}, c_val);
end
disp('Cells of fixed values over nonuniform mesh created successfully!');
for i=1:N_mesh
    c_rand = rand(mesh_uniform{i}.dims);
    c_n{i}= createCellVariable(mesh_nonuniform{i}, c_rand);
end
disp('Cells of random values over nonuniform mesh created successfully!');
%% Part III: create face variables
f_val= 0.5;
% uniform
f_u=cell(N_mesh, 1);
for i=1:N_mesh
    f_u{i}= createFaceVariable(mesh_uniform{i}, f_val);
end
disp('Face variable over uniform mesh created successfully!');
% nonuniform
f_n=cell(N_mesh, 1);
for i=1:N_mesh
    f_n{i}= createFaceVariable(mesh_nonuniform{i}, f_val);
end
disp('Face variable over nonuniform mesh created successfully!');
%% Part IV: Test boundary conditions
BC_n=cell(N_mesh, 1);
for i=1:N_mesh
    BC_n{i}=createBC(mesh_uniform{i});
%% Part IV:
% uniform
figure(1);
for i=1:N_mesh
    subplot(3, 3, i);
    visualizeCells(c_u{i});
end
disp('Cells of random values over uniform mesh visualized successfully!');
% nonuniform
figure(2);
for i=1:N_mesh
    subplot(3, 3, i);
    visualizeCells(c_n{i});
end
disp('Cells of random values over nonuniform mesh visualized successfully!');