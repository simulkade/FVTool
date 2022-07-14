classdef CalculableStructTest < matlab.unittest.TestCase
    properties
        x
        MWg
        x_O2_air
        x_N2_air
    end

    methods(TestMethodSetup)

        function setup(testCase)
            testCase.x_O2_air = 0.2095;
            testCase.x_N2_air = 0.7808;
            testCase.x = CalculableStruct(struct('O2', testCase.x_O2_air, 'N2', testCase.x_N2_air));
        end
    end

    methods 
        function verifySame(testCase, a, b)
            testCase.verifyEqual(a,b, 'RelTol', 1e-6);
        end
    end

    methods(Test)
        function test_construction(testCase)
            testCase.verifySame(testCase.x.S.N2, testCase.x_N2_air);
        end

        function test_properties(testCase)
            testCase.x.O2 = 3;
            testCase.verifySame(testCase.x.O2, 3);
        end
        
        function test_add_field(testCase)
            testCase.x.add_field("AA", 2);
            testCase.verifySame(testCase.x.AA, 2);

            testCase.x.add_field("Ar", 1 - testCase.x.O2 - testCase.x.N2);
            testCase.verifySame(testCase.x.Ar, 1-testCase.x_O2_air-testCase.x_N2_air);
        end

        function test_plus(testCase)
            % Two CalculableStructs
            res = testCase.x + testCase.x;
            testCase.verifySame(res.O2, 2*testCase.x_O2_air);

            % CalculablStruct left and CellTable right
            Nx = 3;
            L = 1;
            m = createMesh1D(Nx, L);
            ct = CellTable.from_struct(m, testCase.x);
            res = testCase.x + ct;
            testCase.verifySame(res.O2.value(2), 2 * testCase.x_O2_air);
            testCase.verifySame(class(res), 'CellTable');

            my_struct = struct();
            my_struct.O2 = 4;
            my_struct.N2 = 2;
            res = my_struct + testCase.x;
            testCase.verifySame(res.O2, testCase.x_O2_air + 4);
            testCase.verifySame(res.N2, testCase.x_N2_air + 2);

            my_struct = struct();
            my_struct.O2 = 4;
            my_struct.N2 = 2;
            res =  testCase.x + my_struct;
            testCase.verifySame(res.O2, testCase.x_O2_air + 4);
            testCase.verifySame(res.N2, testCase.x_N2_air + 2);
        end

        function test_uminus(testCase)
            res = -testCase.x;
            testCase.verifySame(res.O2, -testCase.x_O2_air);
        end

        function test_minus(testCase)
            y.O2 = 0.5;
            y.N2 = 0.5;
            res = y-testCase.x;
            testCase.verifySame(res.O2, 0.5-testCase.x_O2_air);
        end

        function test_times(testCase)
            res = testCase.x .* testCase.x;
            testCase.verifySame(res.O2, testCase.x_O2_air * testCase.x_O2_air);
            
            % CalculablStruct left and CellTable right
            Nx = 3;
            L = 1;
            m = createMesh1D(Nx, L);
            ct = CellTable.from_struct(m, testCase.x);
            res = testCase.x .* ct;
            testCase.verifySame(res.O2.value(2), testCase.x_O2_air * testCase.x_O2_air);
            testCase.verifySame(class(res), 'CellTable');

            my_struct = struct();
            my_struct.O2 = 4;
            my_struct.N2 = 2;
            res = my_struct .* testCase.x;
            testCase.verifySame(res.O2, testCase.x_O2_air .* 4);
            testCase.verifySame(res.N2, testCase.x_N2_air .* 2);

            my_struct = struct();
            my_struct.O2 = 4;
            my_struct.N2 = 2;
            res =  testCase.x .* my_struct;
            testCase.verifySame(res.O2, testCase.x_O2_air .* 4);
            testCase.verifySame(res.N2, testCase.x_N2_air .* 2);
        end

        function test_rdivide(testCase)
            res = testCase.x ./ testCase.x;
            testCase.verifySame(res.O2, 1);
            
            % CalculablStruct left and CellTable right
            Nx = 3;
            L = 1;
            m = createMesh1D(Nx, L);
            ct = CellTable.from_struct(m, testCase.x);
            res = testCase.x ./ ct;
            testCase.verifySame(res.O2.value(2), 1);
            testCase.verifySame(class(res), 'CellTable');

            my_struct = struct();
            my_struct.O2 = 4;
            my_struct.N2 = 2;
            res = my_struct ./ testCase.x;
            testCase.verifySame(res.O2, 4 / testCase.x_O2_air);

            my_struct = struct();
            my_struct.O2 = 4;
            my_struct.N2 = 2;
            res =  testCase.x ./ my_struct;
            testCase.verifySame(res.O2, testCase.x_O2_air / 4);
        end

        function test_copy(testCase)
            other = copy(testCase.x);
            other.O2 = 3;

            testCase.verifySame(testCase.x.O2, testCase.x_O2_air);
        end

        function test_sum(testCase)
            res = sum(testCase.x);
            testCase.verifySame(testCase.x_O2_air + testCase.x_N2_air, res);
        end

        function test_to_vec(testCase)
            vec= testCase.x.toVec();
            testCase.verifySame([testCase.x_N2_air; testCase.x_O2_air], vec);
        end

        function test_patch_vec(testCase)
            vec = [3; 4];
            new_cs = CalculableStruct.patch_vec(testCase.x, vec);
            testCase.verifySame(new_cs.N2, 3);
        end

    end

end

