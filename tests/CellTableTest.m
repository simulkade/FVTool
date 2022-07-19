classdef CellTableTest < matlab.unittest.TestCase
    properties
        x
        m
        x_O2_air
        x_N2_air
    end

    methods(TestMethodSetup)

        function setup(testCase)
            Nx = 3;
            L = 1;
            ms = createMesh1D(Nx, L);
            testCase.m = ms;
            testCase.x_O2_air = 0.2095;
            testCase.x_N2_air = 0.7808;
            O2 = createCellVariable(ms, testCase.x_O2_air);
            N2 = createCellVariable(ms, testCase.x_N2_air);
            testCase.x = CellTable(ms, table(O2, N2));
            testCase.x.equip_prop();
        end
    end

    methods 
        function verifySame(testCase, a, b)
            testCase.verifyEqual(a,b, 'RelTol', 1e-6);
        end
    end

    methods(Test)
        function test_construction(testCase)
            arr = zeros(1, 2, 5);
            arr(1,1,:) = testCase.x_N2_air;
            arr(1,2,:) = testCase.x_O2_air;
            testCase.verifyEqual(testCase.x.A, arr);
        end


        function test_property_get(testCase)
            testCase.verifySame(testCase.x.O2.value(3), testCase.x_O2_air);
        end

        function test_properties_set(testCase)
            testCase.x.O2.value(3) = 3;
            testCase.verifySame(testCase.x.O2.value(3), 3);
        end

        function test_get_cv(testCase)
            testCase.verifySame(testCase.x.get_cv("O2").value(2), testCase.x_O2_air);
        end

        function test_patch_cv(testCase)
            cv_O2 = createCellVariable(testCase.m, 1);
            testCase.x.patch_cv("O2", cv_O2);
            testCase.verifySame(testCase.x.O2.value(2), 1);
        end
        
        function test_add_field(testCase)
            testCase.x.add_field("AA", createCellVariable(testCase.m, 1));
            testCase.x.equip_prop()
            testCase.verifySame(testCase.x.AA.value(2), 1);

            testCase.x.add_field("Ar", 1 - testCase.x.O2 - testCase.x.N2);
            testCase.x.equip_prop();
            testCase.verifySame(testCase.x.Ar.value(2), 1-testCase.x_O2_air-testCase.x_N2_air);
        end

        function from_array(testCase)
            ct = CellTable.from_array(testCase.x.mesh, testCase.x.A, testCase.x.field_struct);
            ct.equip_prop();
            testCase.verifySame(testCase.x.O2.value(2), testCase.x_O2_air);
        end

        function test_from_calculable_struct(testCase)
            Nx = 3;
            L = 1;
            ms = createMesh1D(Nx, L);

            my_struct = struct();
            my_struct.O2 = 4;
            my_struct.N2 = 2;
            cs = CalculableStruct(my_struct);
            
            ct = CellTable.from_calculable_struct(ms, cs);
            ct.equip_prop();
            testCase.verifySame(ct.O2.value(1), 4);

            % my_struct = struct();
            % my_struct.O2 = (1:3)';
            % my_struct.N2 = 2;
            % cs = CalculableStruct(my_struct);
            
            % ct = CellTable.from_struct(ms, my_struct);
            % ct.equip_prop();
            % testCase.verifySame(ct.O2.value(2), 1);
        end

        function test_plus(testCase)
            res = testCase.x + testCase.x;
            res.equip_prop();
            testCase.verifySame(res.O2.value(1), 2*testCase.x_O2_air);

            cell_var = createCellVariable(testCase.m, (1:3)');
            res = testCase.x + cell_var;
            res.equip_prop();
            testCase.verifySame(res.O2.value(1 + 2), testCase.x_O2_air + 2);

            cs = CalculableStruct(struct("O2", 1, "N2", 1));
            res = testCase.x + cs;
            res.equip_prop();
            testCase.verifySame(res.O2.value(1), testCase.x_O2_air + 1);

            cs = CalculableStruct(struct("O2", 1, "N2", 1));
            res = cs + testCase.x;
            res.equip_prop();
            testCase.verifySame(res.O2.value(1), testCase.x_O2_air + 1);

            res = testCase.x + 3;
            res.equip_prop();
            testCase.verifySame(res.O2.value(1), testCase.x_O2_air + 3);

            res = 3 + testCase.x;
            res.equip_prop();
            testCase.verifySame(res.O2.value(1), testCase.x_O2_air + 3);
        end

        function test_uminus(testCase)
            res = -testCase.x;
            res.equip_prop();
            testCase.verifySame(res.O2.value(1), -testCase.x_O2_air);
        end

        function test_minus(testCase)
            res = testCase.x - testCase.x;
            res.equip_prop();
            testCase.verifySame(res.O2.value(1), 0);

            cs = CalculableStruct(struct("O2", 1, "N2", 2));
            res = testCase.x - cs;
            res.equip_prop();
            testCase.verifySame(res.O2.value(1), testCase.x_O2_air - 1);

            cs = CalculableStruct(struct("O2", 1, "N2", 2));
            res = cs - testCase.x;
            res.equip_prop();
            testCase.verifySame(res.O2.value(1), 1 - testCase.x_O2_air);

            res = testCase.x - 3;
            res.equip_prop();
            testCase.verifySame(res.O2.value(1), testCase.x_O2_air - 3);

            res = 3 - testCase.x;
            res.equip_prop();
            testCase.verifySame(res.O2.value(1), 3 - testCase.x_O2_air);
        end

        function test_times(testCase)
            res = testCase.x .* testCase.x;
            res.equip_prop();
            testCase.verifySame(res.O2.value(1), testCase.x_O2_air * testCase.x_O2_air);

            my_struct = struct();
            my_struct.O2 = 4;
            my_struct.N2 = 2;
            cs = CalculableStruct(my_struct);

            res = testCase.x .* cs;
            res.equip_prop();
            testCase.verifySame(res.O2.value(1), testCase.x_O2_air * 4);

            res = cs .* testCase.x;
            res.equip_prop();
            testCase.verifySame(res.O2.value(1), testCase.x_O2_air * 4);

            res = testCase.x .* 3;
            res.equip_prop();
            testCase.verifySame(res.O2.value(1), testCase.x_O2_air * 3);

            res = 3 .* testCase.x;
            res.equip_prop();
            testCase.verifySame(res.O2.value(1), testCase.x_O2_air * 3);
        end

        function test_rdivide(testCase)
            res = testCase.x ./ testCase.x;
            res.equip_prop();
            testCase.verifySame(res.O2.value(1), 1);

            my_struct = struct();
            my_struct.O2 = 4;
            my_struct.N2 = 2;
            cs = CalculableStruct(my_struct);

            res = testCase.x ./ cs;
            res.equip_prop();
            testCase.verifySame(res.O2.value(1), testCase.x_O2_air / 4);

            res = cs ./ testCase.x;
            res.equip_prop();
            testCase.verifySame(res.O2.value(1), 4 / testCase.x_O2_air);

            res = 2 ./ testCase.x;
            res.equip_prop();
            testCase.verifySame(res.O2.value(1), 2 / testCase.x_O2_air);

            res = testCase.x ./ 2;
            res.equip_prop();
            testCase.verifySame(res.O2.value(1), testCase.x_O2_air / 2);

        end

        function test_power(testCase)
            res = testCase.x.^2;
            res.equip_prop();
            testCase.verifySame(res.O2.value(2), testCase.x_O2_air^2); 
        end


        function test_sum(testCase)
            testCase.x.add_field("Ar", 1 - testCase.x.O2 - testCase.x.N2);

            sum_x = sum(testCase.x);
            testCase.verifySame(sum_x.value(2), 1);
        end

        function test_mass_mol(testCase)
            testCase.x.add_field("Ar", 1 - testCase.x.O2 - testCase.x.N2);

            MWg = CalculableStruct(struct('N2', 28e-3, 'O2', 32e-3, 'Ar', 40e-3));
            MWg.equip_prop();
            x_Ar_air = 1-testCase.x_O2_air - testCase.x_N2_air;
            manual_cal = (testCase.x_O2_air * MWg.O2) / (testCase.x_O2_air * MWg.O2 + ...
                testCase.x_N2_air * MWg.N2 + x_Ar_air * MWg.Ar);

            y = testCase.x .* MWg ./ (sum(testCase.x .* MWg));
            y.equip_prop();
            testCase.verifySame(y.O2.value(2), manual_cal);
            testCase.verifyEqual(y.O2.value(2), 0.23, 'RelTol', 1e-2);
        end

    end

end

