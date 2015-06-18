function [maxprop,uncor_prop,maxcsize] = meta_stochastic_activation(MC_Setup)
% MC_Setup = meta_stochastic_activation(DB)
%
% This function randomizes n locations for each study within analysis mask
% and computes null-hypothesis weighted proportion of activated studies
% (contrasts)
% for whole-brain FWE corrected Monte Carlo
% 
% This function is part of Meta_Activation_Fwe
%
% Requesting uncor_prop and maxcsize slows processing.
%
% Tor Wager, May 2006


nc = length(MC_Setup.n);
v = size(MC_Setup.xyzlist,1);

ivectors = zeros(v,1);       % initialize output

str1 = sprintf('Convolving: %03d',0);,fprintf(1,str1);

for c = 1:nc

    fprintf(1,'\b\b\b%03d',c);
    
    % get convolved, randomized coords for this study
    indic = meta_stochastic_mask(MC_Setup.xyzlist,MC_Setup.n(c),MC_Setup.r);

    % weight and add to summed image vector
    indic = indic .* MC_Setup.wts(c);
    ivectors = ivectors + indic;
    
end

erase_string(str1);

% ------------------------------------------------
% save summary info from this map only
% ------------------------------------------------
str1 = sprintf('Saving summary info');,fprintf(1,str1);

maxprop = max(ivectors);

if nargout > 1
    uncor_prop = prctile(ivectors,[95 99 99.9]);

    if nargout > 2
        for i = 1:length(uncor_prop)

            maxcsize(i) = get_max_cluster_size(ivectors, MC_Setup.xyzlist,uncor_prop(i));

        end
    end
end

erase_string(str1);


return


function  maxcsize = get_max_cluster_size(ivectors,xyzlist,prop)
% save max cluster size info
wh = find(ivectors >= prop);

if length(wh) > 50000, maxcsize = Inf;, end

xyz = xyzlist(wh,:)';
clust = spm_clusters(xyz);
u = unique(clust);
for i = 1:length(u)
    csize(i) = sum(clust == i);
end
maxcsize = max(csize);
return



function erase_string(str1)
fprintf(1,repmat('\b',1,length(str1))); % erase string
%fprintf(1,'\n');
return
