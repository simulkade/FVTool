classdef MeshStructure
    %CELLVARIABLE Summary of this class goes here
    %   Detailed explanation goes here

    properties
        dimension
        dims
        cellsize
        cellcenters
        facecenters
        corners
        edges
    end

    methods
        function meshVar = MeshStructure(dimension, dims, cellsize, ...
          cellcenters, facecenters, corners, edges)
            if nargin>0
                meshVar.dimension = dimension;
                meshVar.dims = dims;
                meshVar.cellsize = cellsize;
                meshVar.cellcenters = cellcenters;
                meshVar.facecenters = facecenters;
                meshVar.corners= corners;
                meshVar.edges= edges;
            end
        end
    end
end
