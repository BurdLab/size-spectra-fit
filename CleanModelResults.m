function [poc_prefac, poc_slopes, n_prefac, n_slopes] = ...
    CleanModelResults(poc_prefac, poc_slopes, n_prefac, n_slopes)

% Find indices of results with positive prefactors and slopes

indx_pos = find(poc_prefac > 0 & poc_slopes > 0 & n_prefac > 0 & n_slopes > 0); 

% pos_poc_prefac = find(poc_prefac > 0);
% pos_n_prefac   = find(n_prefac > 0);
% pos_prefac     = union(pos_poc_prefac, pos_n_prefac);
% 
% pos_poc_slopes = find(poc_slopes > 0);
% pos_n_slopes   = find(n_slopes > 0);
% pos_slopes     = union(pos_poc_slopes, pos_n_slopes);
% 
% indx_pos = union(pos_prefac, pos_slopes);

% Now filter data

poc_prefac = poc_prefac(indx_pos);
poc_slopes = poc_slopes(indx_pos);
n_prefac   = n_prefac(indx_pos);
n_slopes   = n_slopes(indx_pos);





