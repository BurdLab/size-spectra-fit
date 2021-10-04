function [f, J] = TestFunc2(x, data, p, size_ranges)
% 
% TestFunc1 

a = x(1);
b = x(2);

f = zeros(2,1);

f(1) = p(1)*a/(1-b+p(2)) * (size_ranges(2)^(1-b+p(2)) - size_ranges(1)^(1-b+p(2))) - data(1); 
f(2) = p(1)*a/(1-b+p(2)) * (size_ranges(4)^(1-b+p(2)) - size_ranges(3)^(1-b+p(2))) - data(2);


% Calculate the Jacobian
if nargout > 1 
    J(1,1) = p(1)/(1-b+p(2)) * (size_ranges(2)^(1-b+p(2)) - size_ranges(1)^(1-b+p(2)));
    J(2,1) = p(1)*a/(1-b+p(2))^2 * (size_ranges(2)^(1-b+p(2)) - size_ranges(1)^(1-b+p(2))) + ...
        p(1)*a/(1-b+p(2))*(log(size_ranges(1))*size_ranges(1)^(1-b+p(2)) - ...
            log(size_ranges(2))*size_ranges(2)^(1-b+p(2))); 
    J(2,1) = p(1)/(1-b+p(2)) * (size_ranges(4)^(1-b+p(2)) - size_ranges(3)^(1-b+p(2)));
    J(2,2) = p(1)*a/(1-b+p(2))^2 * (size_ranges(4)^(1-b+p(2)) - size_ranges(3)^(1-b+p(2))) + ...
        p(1)*a/(1-b+p(2))*(log(size_ranges(3))*size_ranges(3)^(1-b+p(2)) - ...
            log(size_ranges(4))*size_ranges(4)^(1-b+p(2))); 
end

end