classdef CellVariable
    %CELLVARIABLE Summary of this class goes here
    %   Detailed explanation goes here
    % Copyright (c) 2012-2016 Ali Akbar Eftekhari
    % See the license file

    properties
        domain
        value
    end

    methods
        function cv = CellVariable(meshVar, cellval)
            if nargin>0
                cv.domain = meshVar;
                cv.value = cellval; % does not do any dim check!
            end
        end
    end
end
