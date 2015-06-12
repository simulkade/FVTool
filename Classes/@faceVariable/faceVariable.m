classdef faceVariable
    properties
        dimension
        meshsize
        xvalue
        yvalue
        zvalue
        cellsize
    end
    
    methods
        function fv = faceVariable(meshVar, faceval)
            if nargin>0
                fv.dimension = meshVar.dimension;
                fv.meshsize = meshVar.numberofcells;
                fv.cellsize = meshVar.cellsize;
                if nargin==1
                    faceval = ones(floor(meshVar.dimension),1);
                end
                mn = meshVar.numberofcells;    
                if (meshVar.dimension ==1) || (meshVar.dimension==1.5)
                    fv.xvalue = faceval.*ones(mn+1, 1);
                elseif (meshVar.dimension == 2) || (meshVar.dimension == 2.5)
                    fv.xvalue = faceval(1).*ones(mn(1)+1, mn(2));
                    fv.yvalue = faceval(2).*ones(mn(1), mn(2)+1);
                elseif meshVar.dimension == 3
                    fv.xvalue = faceval(1).*ones(mn(1)+1, mn(2), mn(3));
                    fv.yvalue = faceval(2).*ones(mn(1), mn(2)+1, mn(3));
                    fv.zvalue = faceval(3).*ones(mn(1), mn(2), mn(3)+1);
                end
            end
        end
    end
    
end

