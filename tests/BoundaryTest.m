classdef BoundaryTest < matlab.unittest.TestCase

    methods(Test)

        function test_robin(testCase)
            y_fumes.O2 = 0.1;
            y_fumes.N2 = 0.8;
            rob = robin(0, 1, y_fumes); 
            testCase.verifyEqual(rob.c.O2, 0.1);
            testCase.verifyEqual(rob.b, 1);
        end
    end

end

