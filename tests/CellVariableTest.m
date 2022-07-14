classdef CellVariableTest < matlab.unittest.TestCase
    properties
        x
        m
        cvar
    end

    methods(TestMethodSetup)

        function setup(testCase)
            Nx = 4;
            L = 1;
            ms = createMesh1D(Nx, L);
            testCase.m = ms;
            testCase.cvar = createCellVariable(ms, (1:4)');
        end
    end

    methods 
        function verifySame(testCase, a, b)
            testCase.verifyEqual(a,b, 'RelTol', 1e-6);
        end
    end

    methods(Test)
        function test_ival(testCase)
            testCase.verifySame(testCase.cvar.ival(2), 2);
            testCase.cvar.ival = (4:7)';
            testCase.verifySame(testCase.cvar.ival(2), 5);
        end

        function test_applyBC(testCase)
            BC = createBC(testCase.m);
            BC.left = robin(1, 0, 1);
            BC.right = robin(0, 1, 1);
            testCase.cvar = applyBC(testCase.cvar, BC);
            grad = gradientTerm(testCase.cvar); 
            testCase.verifySame(grad.xvalue(1), 1);
            testCase.verifySame(testCase.cvar.right, 1);
        end

    end

end

