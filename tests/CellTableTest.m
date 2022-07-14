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
            testCase.x = CellTable(table(O2, N2));
        end
    end

    methods 
        function verifySame(testCase, a, b)
            testCase.verifyEqual(a,b, 'RelTol', 1e-6);
        end
    end

    methods(Test)
        function test_construction(testCase)
            testCase.verifySame(testCase.x.T.N2.ival(1), testCase.x_N2_air);
        end

        function test_properties(testCase)
            testCase.x.O2.value(3) = 3;
            testCase.verifySame(testCase.x.O2.value(3), 3);
        end
        
        function test_add_field(testCase)
            testCase.x.add_field("AA", createCellVariable(testCase.m, 1));
            testCase.verifySame(testCase.x.AA.value(2), 1);

            testCase.x.add_field("Ar", 1 - testCase.x.O2 - testCase.x.N2);
            testCase.verifySame(testCase.x.Ar.value(2), 1-testCase.x_O2_air-testCase.x_N2_air);
        end

        function test_plus(testCase)
            res = testCase.x + testCase.x;
            testCase.verifySame(res.O2.value(1), 2*testCase.x_O2_air);

            cell_var = createCellVariable(testCase.m, (1:3)');
            res = testCase.x + cell_var;
            testCase.verifySame(res.O2.value(1 + 2), testCase.x_O2_air + 2);

            my_struct = struct();
            my_struct.O2 = 4;
            my_struct.N2 = 2;
            res = my_struct + testCase.x;
            testCase.verifySame(res.O2.value(1), testCase.x_O2_air + 4);
            testCase.verifySame(res.N2.value(1), testCase.x_N2_air + 2);

            my_struct = struct();
            my_struct.O2 = 4;
            my_struct.N2 = 2;
            res =  testCase.x + my_struct;
            testCase.verifySame(res.O2.value(1), testCase.x_O2_air + 4);
            testCase.verifySame(res.N2.value(1), testCase.x_N2_air + 2);

            res = 5 + testCase.x;
            testCase.verifySame(res.O2.value(1), testCase.x_O2_air + 5);
            testCase.verifySame(res.N2.value(1), testCase.x_N2_air + 5);

            res = testCase.x + 5;
            testCase.verifySame(res.O2.value(1), testCase.x_O2_air + 5);
            testCase.verifySame(res.N2.value(1), testCase.x_N2_air + 5);
        end

        function test_uminus(testCase)
            res = -testCase.x;
            testCase.verifySame(res.O2.value(1), -testCase.x_O2_air);
        end

        function test_minus(testCase)
            res = testCase.x - testCase.x;
            testCase.verifySame(res.O2.value(1), 0);

            res = testCase.x - 3;
            testCase.verifySame(res.O2.value(1), testCase.x_O2_air - 3);
        end

        function test_times(testCase)
            res = testCase.x .* testCase.x;
            testCase.verifySame(res.O2.value(1), testCase.x_O2_air * testCase.x_O2_air);

            my_struct = struct();
            my_struct.O2 = 4;
            my_struct.N2 = 2;

            res = testCase.x .* my_struct;
            testCase.verifySame(res.O2.value(1), testCase.x_O2_air * 4);

            res = my_struct .* testCase.x;
            testCase.verifySame(res.O2.value(1), testCase.x_O2_air * 4);

            res = testCase.x .* 3;
            testCase.verifySame(res.O2.value(1), testCase.x_O2_air * 3);

            res = 3 .* testCase.x;
            testCase.verifySame(res.O2.value(1), testCase.x_O2_air * 3);
        end

        function test_rdivide(testCase)
            res = testCase.x ./ testCase.x;
            testCase.verifySame(res.O2.value(1), 1);

            my_struct = struct();
            my_struct.O2 = 4;
            my_struct.N2 = 2;

            res = testCase.x ./ my_struct;
            testCase.verifySame(res.O2.value(1), testCase.x_O2_air / 4);

        end

        function test_from_struct(testCase)
            Nx = 3;
            L = 1;
            ms = createMesh1D(Nx, L);

            my_struct = struct();
            my_struct.O2 = 4;
            my_struct.N2 = 2;
            
            ct = CellTable.from_struct(ms, my_struct);
            testCase.verifySame(ct.O2.value(1), 4);

            my_struct = struct();
            my_struct.O2 = (1:3)';
            my_struct.N2 = 2;
            
            ct = CellTable.from_struct(ms, my_struct);
            testCase.verifySame(ct.O2.value(2), 1);
        end

        function test_sum(testCase)
            testCase.x.add_field("Ar", 1 - testCase.x.O2 - testCase.x.N2);

            sum_x = sum(testCase.x);
            testCase.verifySame(sum_x.value(2), 1);
        end

        function test_mass_mol(testCase)
            testCase.x.add_field("Ar", 1 - testCase.x.O2 - testCase.x.N2);

            MWg = struct('N2', 28e-3, 'O2', 32e-3, 'Ar', 40e-3);
            x_Ar_air = 1-testCase.x_O2_air - testCase.x_N2_air;
            manual_cal = (testCase.x_O2_air * MWg.O2) / (testCase.x_O2_air * MWg.O2 + ...
                testCase.x_N2_air * MWg.N2 + x_Ar_air * MWg.Ar);

            y = testCase.x .* MWg ./ (sum(testCase.x .* MWg));
            testCase.verifySame(y.O2.value(2), manual_cal);
            testCase.verifyEqual(y.O2.value(2), 0.23, 'RelTol', 1e-2);
        end

        function test_apply_BC(testCase)
            bc = createBC(testCase.m);
            x_out.O2 = 0;
            x_out.N2 = 1;
            bc.right = robin(0, 1, x_out);
            
            testCase.x = testCase.x.apply_BC(bc);
            testCase.verifySame(testCase.x.O2.value(end), ...
                testCase.x.O2.ival(end) + 2 * (x_out.O2 - testCase.x.O2.ival(end)));

        end

        function test_fieldsCompatible(testCase)
            a.O2 = 3;
            b.O2 = 2;
            b.AZE = 'okje';

            testCase.verifyError(@() CellTable.fieldsCompatible(a,b),...
                "CellTable:FieldsNotCompatible")
        end

        function test_toArray(testCase)
            mat = testCase.x.toArray();
            testCase.verifySame(mat(1,1), testCase.x_N2_air);
        end
            
    end

end

