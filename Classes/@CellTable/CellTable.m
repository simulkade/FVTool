classdef (InferiorClasses = {?CellVariable, ?CalculableStruct}) CellTable < dynamicprops
    % WIP: the idea is to get a data strcture to hold several labeled
    % CellVariables that one may want to treat in the same way and make some
    % algebraic calculations on. Primary use case is concentration of several
    % species.
    % It is kind of a table but not as efficient, but adapated to manipulate CellVariables

    properties
        A  % inner array
        field_struct  % ordered list
        mesh
    end

    properties (Dependent)
        fields
        nxs
        nf
        T
    end

    methods
        function ct = CellTable(mesh, inner_table)
            arguments
                mesh MeshStructure
                inner_table table = table()
            end
            % ct.T = inner_table(:, sort(inner_table.Properties.VariableNames));
            fields = string(sort(inner_table.Properties.VariableNames))';
            ct.mesh = mesh;
            Nxs = ct.mesh.dims(1) + 2;
            ct.A = zeros(1, numel(fields), Nxs); 
            ct.field_struct = struct();

            for idx = 1:numel(fields)
                field = fields(idx);
                ct.field_struct.(field) = idx;
                ct.A(1, idx, :) = reshape(inner_table.(field).value, [1 1 Nxs]); 
            end

        end

        function add_prop(self, name)
            if ~isprop(self, name)
                prop = addprop(self, name);
                prop.Dependent = true;
                prop.GetMethod = @(obj)CellTable.getDynamicProp(obj,name);
                prop.SetMethod = @(obj,val)CellTable.setDynamicProp(obj,name,val);
            end
        end

        function equip_prop(self)
            for idx = 1:numel(self.fields)
                field = self.fields(idx);
                self.add_prop(field);
            end
        end

        function add_field(self, name, cell_variable)
            old_arr = self.A;
            old_field_struct = self.field_struct;
            new_fields = sort([self.fields; name]);
            self.field_struct = struct();
            Nxs = self.mesh.dims(1) + 2;
            self.A = zeros(1, numel(new_fields), Nxs); 
            for idx = 1:numel(new_fields)
                field = new_fields(idx);
                self.field_struct.(field) = idx;
                if field == name
                    self.A(1, idx, :) = reshape(cell_variable.value, [1 1 Nxs]); 
                else
                    self.A(1, idx, :) = old_arr(1, old_field_struct.(field), :);
                end
            end
        end

        function f = get.fields(self)
            f = string(fieldnames(self.field_struct));
        end

        function tab = get.T(self)
            tab = table();
            for idx = 1:self.nf
                field = self.fields(idx);
                tab.(field) = reshape(self.A(1, idx, :), [self.nxs 1]);
            end
        end

        function repr(self)
            disp(self.T);
        end

        function res = sum(self)
            % Does horizontal sum
            v = reshape(sum(self.A, 2), [self.mesh.dims(1)+2 1]);
            res = createCellVariable(self.mesh, v(2:end-1));
        end

        function cv = get_cv(self, name)
            vec = reshape(self.A(1, self.field_struct.(name), :), [self.nxs 1]);
            cv = createCellVariable(self.mesh, vec(2:end-1));
        end

        function patch_cv(self, name, cv)
            self.A(1, self.field_struct.(name), :) = reshape(cv.value, [1 1 self.nxs]);
        end

        function patch_vec(self, name, vec)
            self.A(1, self.field_struct.(name), :) = reshape(vec, [1 1 self.nxs]);
        end

        function patch_ct(self, ct)
            self.A = ct.A;
        end

        function n = get.nxs(self)
            n = self.mesh.dims(1) + 2;
        end

        function n = get.nf(self)
            n = numel(self.fields);
        end

        function new_obj = copy(obj)
            new_obj = CellTable.from_array(obj.mesh, obj.A, obj.field_struct);
        end

        function obj = toArray(self)
            len = numel(self.fields);
            obj = zeros(self.mesh.dims(1) + 2, len);

            for idx_f = 1:len
                field = self.fields{idx_f};
                obj(:,idx_f) = self.(field).value;
            end
        end
            
            
    end

    methods (Static)

        function ct = from_array(mesh, arr, field_struct)
            ct = CellTable(mesh);
            ct.A = arr;
            ct.field_struct = field_struct;
        end


        function f = safe_fields(p)
            if class(p) == "struct"
                f = sort(string(fieldnames(p)))';
            else
                f = p.fields;
            end
        end

        function fieldsCompatible(p,q)
            a = CellTable.safe_fields(p);
            b = CellTable.safe_fields(q);
            if (numel(setdiff(a,b)) > 0 || numel(setdiff(b,a)) > 0)
                throw(MException("CellTable:FieldsNotCompatible", "Operations need fields to be the same"));
            end
        end

        function val = getDynamicProp(obj,name)
          idx = obj.field_struct.(name);
          val = createCellVariable(obj.mesh, reshape(obj.A(1, idx, 2:end-1), [obj.mesh.dims(1) 1]));
        end

        function setDynamicProp(obj,name,val)
            idx = obj.field_struct.(name);
            if class(val) == "CellVariable"
                Nxs = val.domain.dims(1) + 2;
                obj.A(1,idx,:) = reshape(val.value, [1 1 Nxs]);
            else 
                Nxs = numel(val);
                obj.A(1,idx,:) = reshape(val, [1 1 Nxs]);
            end
        end


        function ct = from_calculable_struct(mesh, cs)
            ct = CellTable(mesh);
                for idx = 1:numel(cs.fields) 
                    ct.add_field(cs.fields(idx), createCellVariable(mesh, cs.V(idx)));
                end
        end

    end

end

