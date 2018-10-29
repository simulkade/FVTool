classdef BoundaryCondition
    %CELLVARIABLE Summary of this class goes here
    %   Detailed explanation goes here
    % Copyright (c) 2012-2016 Ali Akbar Eftekhari
    % See the license file

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
