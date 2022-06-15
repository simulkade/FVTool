classdef CellVariable
    %CellVariable class 
    % 

    properties
        domain
        value
    end

    properties (Dependent)
        ival
        left
        right
    end



    methods
        function cv = CellVariable(meshVar, cellval)
            if nargin>0
                cv.domain = meshVar;
                cv.value = cellval; % does not do any dim check!
            end
        end

        function self = apply_BC(self, BC)
            self.value = cellBoundary(self.ival, BC);
        end

        function r = get.ival(self)
            % Inner value
            r = self.value(2:end-1);
        end

        function r = get.left(self)
            % Value at left boundary
            r = (self.value(1) + self.value(2))/2;
        end

        function r = get.right(self)
            % Value at left boundary
            r = (self.value(end) + self.value(end-1))/2;
        end
    end
end
