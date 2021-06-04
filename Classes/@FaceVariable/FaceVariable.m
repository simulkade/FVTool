classdef FaceVariable
% FaceVariable class
    properties
        domain
        xvalue
        yvalue
        zvalue
    end

    methods
        function fv = FaceVariable(meshVar, facevalX, facevalY, facevalZ)
            if nargin>0
                fv.domain = meshVar;
                fv.xvalue = facevalX;
                fv.yvalue = facevalY;
                fv.zvalue = facevalZ;
            end
        end
    end

end
