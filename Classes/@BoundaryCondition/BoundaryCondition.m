classdef BoundaryCondition
    %CELLVARIABLE Summary of this class goes here
    %   Detailed explanation goes here

    properties
        domain
        left
        right
        bottom
        top
        back
        front
    end

    methods
        function BC = BoundaryCondition(meshvar, left, right, bottom, top, back, front)
            if nargin>0
                BC.domain=meshvar;
                BC.left=left;
                BC.right=right;
                BC.bottom=bottom;
                BC.top=top;
                BC.back=back;
                BC.front=front;
            end
        end
    end
end
