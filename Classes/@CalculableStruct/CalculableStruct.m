classdef CalculableStruct < dynamicprops
    % WIP: Looks like a lot of work (and probably bad performance) 
    % just to allow nice labeled calculation but trying
    properties
        S  % inner table
    end

    properties (Dependent)
        fields
    end

    methods
        function cs = CalculableStruct(inner_struct)
            arguments
                inner_struct struct = struct()
            end
            inner_struct = orderfields(inner_struct);
            cs.S = inner_struct;
            fields = fieldnames(inner_struct);

            for idx = 1:numel(fields)
                field = fields{idx};
                cs.add_prop_from_struct(field);
            end
        end

        function add_prop_from_struct(self, name)
                prop = addprop(self, name);
                prop.Dependent = true;
                prop.GetMethod = @(obj)CalculableStruct.getDynamicProp(obj,name);
                prop.SetMethod = @(obj,val)CalculableStruct.setDynamicProp(obj,name,val);
        end

        function add_field(self, name, value)
            self.S.(name) = value;
            self.add_prop_from_struct(name);
        end

        function f = get.fields(self)
            f = fieldnames(self.S);
        end

        function str = repr(self)
            str = disp(self.S);
        end

        function res = sum(self)
            % Does horizontal sum
            res = 0;
            for idx = 1:numel(self.fields)
                field = self.fields{idx};
                res = res + self.(field);
            end
        end

        function new_obj = copy(self)
            new_obj = CalculableStruct(self.S);
        end

        function vec = toVec(self)
            len = numel(self.fields);
            vec = zeros(len,1);
            for idx_f = 1:len
                field = self.fields{idx_f};
                vec(idx_f) = self.(field);
            end
        end


    end

    methods (Static)

        function val = getDynamicProp(obj,name)
          val = obj.S.(name);
        end

        function setDynamicProp(obj,name,val)
          obj.S.(name) = val;
        end

        function obj = patch_vec(cs, vec)
            obj = copy(cs);
            for idx_f = 1:numel(cs.fields)
                field = cs.fields{idx_f};
                obj.S.(field) = vec(idx_f);
            end
        end

    end

end

