classdef StructUtilitiesTest < matlab.unittest.TestCase
    properties
    end

    methods(TestMethodSetup)
    end

    methods 
    end

    methods(Test)
        function structProd_test(testCase)
            a.O2 = 1;
            a.N2 = 2;
            b.O2 = 1;
            b.N2 = 2;
            c = structProd([a b]);
            testCase.verifyEqual(c.O2, 1);
            testCase.verifyEqual(c.N2, 4);
        end

        function structScalarProd(testCase)
            a.O2 = 1;
            a.N2 = 2;
            b = 4;
            c = structScalarProd(a, b);
            testCase.verifyEqual(c.O2, 4);
            testCase.verifyEqual(c.N2, 8);
        end

        function structSum(testCase)
            a.O2 = 1;
            a.N2 = 2;
            c = structSum(a);
            testCase.verifyEqual(c, 3);
        end
    end
end
