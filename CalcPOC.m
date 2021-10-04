function poc = CalcPOC(poc_slopes, poc_prefac, size_ranges)

poc = zeros(length(poc_slopes), 2);

% Calculate objective function
poc(:,1) = poc_prefac./(1-poc_slopes) .* (size_ranges(2).^(1-poc_slopes) - size_ranges(1).^(1-poc_slopes)); 
poc(:,2) = poc_prefac./(1-poc_slopes) .* (size_ranges(4).^(1-poc_slopes) - size_ranges(3).^(1-poc_slopes));
