function [f, J] = TestFunc1(x, data, size_ranges)
% 
% TestFunc1 

a = x(1);
b = x(2);

f = zeros(2,1);

% Calculate objective function
f(1) = a/(1-b) * (size_ranges(2)^(1-b) - size_ranges(1)^(1-b)) - data(1); 
f(2) = a/(1-b) * (size_ranges(4)^(1-b) - size_ranges(3)^(1-b)) - data(2);

% Calculate Jacobian
if nargout > 1
    J(1,1) = 1/(1-b)*(size_ranges(2)^(1-b) - size_ranges(1)^(1-b));
    J(1,2) = a/(1-b)^2 * (size_ranges(2)^(1-b) - size_ranges(1)^(1-b)) + ...
        a/(1-b) * (log(size_ranges(1))*size_ranges(1)^(1-b) - ...
        log(size_ranges(2))*size_ranges(2)^(1-b));
    J(2,1) = 1/(1-b)*(size_ranges(4)^(1-b) - size_ranges(3)^(1-b));
    J(2,2) = a/(1-b)^2 * (size_ranges(4)^(1-b) - size_ranges(3)^(1-b)) + ...
        a/(1-b) * (log(size_ranges(3))*size_ranges(3)^(1-b) - ...
        log(size_ranges(4))*size_ranges(4)^(1-b));
end

end