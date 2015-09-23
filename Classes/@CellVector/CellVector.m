classdef CellVector
    properties
        domain
        xvalue
        yvalue
        zvalue
    end

    methods
        function fv = CellVector(meshVar, facevalX, facevalY, facevalZ)
            if nargin>0
                fv.domain = meshVar;
                fv.xvalue = facevalX;
                fv.yvalue = facevalY;
                fv.zvalue = facevalZ;
            end
        end
    end

end
