classdef testclass2 < testclass1
    methods
        function o = testclass2(a)
            o@testclass1(a);
            o.b = 1;
        end
    end
end