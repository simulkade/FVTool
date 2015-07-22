classdef CellVariable
    %CELLVARIABLE Summary of this class goes here
    %   Detailed explanation goes here

    properties
        domain
        value
    end

    methods
        function cv = CellVariable(meshVar, cellval)
            if nargin>0
                cv.domain = meshVar;
                cv.value = cellval.*ones(meshVar.dims);
            end
        end
    end
end
