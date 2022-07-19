classdef (InferiorClasses = {?CellVariable}) CalculableStruct < dynamicprops
    % WIP: Looks like a lot of work (and probably bad performance) 
    % just to allow nice labeled calculation but trying
    properties
        V  % inner vector
        field_struct  % records the field addresses
    end

    properties (Dependent)
        fields
        struct_repr
        nf
    end

    methods
        function cs = CalculableStruct(base_struct)
            arguments
                base_struct struct = struct()
            end
            base_struct = orderfields(base_struct);
            fields = fieldnames(base_struct);
            N = numel(fields);
            cs.V = zeros(N, 1);
            for idx = 1:numel(fields)
                field = fields{idx};
                cs.V(idx) = base_struct.(field);
                cs.field_struct.(field) = idx;
            end
        end

        function equip_prop(self)
            for idx = 1:numel(self.V)
                field = self.fields(idx);
                self.add_prop(field);
            end
        end

        function add_prop(self, name)
            if ~isprop(self, name)
                prop = addprop(self, name);
                prop.Dependent = true;
                prop.GetMethod = @(obj)CalculableStruct.getDynamicProp(obj,name);
                prop.SetMethod = @(obj,val)CalculableStruct.setDynamicProp(obj,name,val);
            end
        end

        function self = add_field(self, name, value)
            st = self.struct_repr;
            st.(name) = value;
            base_struct = orderfields(st);
            new_fields = fieldnames(base_struct);
            N = numel(new_fields);
            self.V = zeros(N, 1);
            self.field_struct = struct();
            for idx = 1:numel(new_fields)
                field = new_fields{idx};
                self.V(idx) = base_struct.(field);
                self.field_struct.(field) = idx;
            end
            self.add_prop(name);
        end

        function f = get.fields(self)
            f = string(fieldnames(self.field_struct));
        end

        function n = get.nf(self)
            n = numel(self.fields);
        end

        function st = get.struct_repr(self)
            st = struct();
            for idx = 1:numel(self.fields)
                field = self.fields(idx);
                st.(field) = self.V(idx);
            end
        end

        function str = repr(self)
            str = disp(self.struct_repr);
        end

        function res = sum(self)
            res = sum(self.V);
        end

        function new_obj = copy(self)
            new_obj = CalculableStruct.from_vec(self.V, self.field_struct);
        end
    end

    methods (Static)

        function obj = from_vec(vec, field_s)
            obj = CalculableStruct();
            obj.V = vec;
            obj.field_struct = field_s;
        end

        function val = getDynamicProp(obj,name)
          val = obj.V(obj.field_struct.(name));
        end

        function setDynamicProp(obj,name,val)
          obj.V(obj.field_struct.(name)) = val;
        end

    end

end

