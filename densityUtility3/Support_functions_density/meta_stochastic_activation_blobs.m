function [maxprop,uncor_prop,maxcsize] = meta_stochastic_activation_blobs(MC_Setup)
% MC_Setup = meta_stochastic_activation_blobs(DB)
%
% This function randomizes n contiguous blobs for each study within analysis mask
% and computes null-hypothesis weighted proportion of activated studies
% (contrasts)
% for whole-brain FWE corrected Monte Carlo
%
% If we have ctxtxi field describing contrasts (set up in meta_prob_activation),
% save summary stats by condition as well
%
% This function is part of Meta_Activation_FWE
%
% Requesting uncor_prop and maxcsize slows processing.
%
% Tor Wager, May 2006

% Checks
if any(isnan(MC_Setup.wts))
    error('MC_Setup.wts has NaNs. Fix before running');
end

% If we have ctxtxi field describing contrasts, save summary stats by
% condition as well
do_by_condition = isfield(MC_Setup,'ctxtxi');

nc = length(MC_Setup.n);
v = size(MC_Setup.volInfo.xyzlist,1);

ivectors = zeros(v,1);       % initialize output

if do_by_condition, ivectors2 = sparse(v,nc); end  % for contrasts

str1 = sprintf('Shuffling blob centers: %03d',0); fprintf(1,str1);

for c = 1:nc
    
    if mod(c,10) == 0, fprintf(1,'\b\b\b%03d',c); end
    
    % Make image with random blob locations for each study
    % ------------------------------------------------
    indic = meta_stochastic_mask_blobs(MC_Setup.volInfo,MC_Setup.cl{c});
    
    if do_by_condition, ivectors2(:,c) = indic; end
    
    % weight and add to summed image vector (for speed)
    indic = indic .* MC_Setup.wts(c);
    ivectors = ivectors + indic;
    
    
end

% Get contrast values
% ------------------------------------------------
if do_by_condition
    % contrasts x voxels
    contrast_est = (MC_Setup.ctxtxi * ivectors2')';
end

erase_string(str1);

% ------------------------------------------------
% save summary info from this map only
% ------------------------------------------------
str1 = sprintf('Saving summary info'); fprintf(1,str1);

maxprop.act = max(ivectors);
if do_by_condition
    maxprop.poscon = max(contrast_est);
    maxprop.negcon = min(contrast_est);
end

if nargout > 1
    uncor_prop.act = prctile(ivectors,[95 99 99.9]);
    
    if do_by_condition
        uncor_prop.poscon = prctile(contrast_est,[95 99 99.9]);
        uncor_prop.negcon = prctile(contrast_est,[5 1 .1]);
        
        ncons = size(contrast_est,2);
        % transpose if multiple contrasts entered -- prctile returns col
        % vectors for multiple columns
        if ncons > 1
            uncor_prop.poscon = uncor_prop.poscon';
            uncor_prop.negcon = uncor_prop.negcon';
        end
    end
    
    if nargout > 2
        for i = 1:length(uncor_prop.act)
            
            maxcsize.act(i) = get_max_cluster_size(ivectors, MC_Setup.volInfo.xyzlist,uncor_prop.act(i));
            
            if do_by_condition
                % uncor.prop.poscon : rows are contrasts, cols are
                % different thresholds
                
                for j = 1:ncons
                    maxcsize.poscon(j,i) = get_max_cluster_size(contrast_est(:,j), MC_Setup.volInfo.xyzlist,uncor_prop.poscon(j,i));
                    maxcsize.negcon(j,i) = get_max_cluster_size(-contrast_est(:,j), MC_Setup.volInfo.xyzlist,-uncor_prop.negcon(j,i));
                end
            end
        end
    end
end

erase_string(str1);


return


function  maxcsize = get_max_cluster_size(ivectors,xyzlist,prop)

% save max cluster size info
wh = find(ivectors >= prop);

if length(wh) > 50000, maxcsize = Inf; end

xyz = xyzlist(wh,:)';
clust = spm_clusters(xyz);
u = unique(clust);
clsize = [];         % init to avoid error in Matlab 2013/csize is a function
for i = 1:length(u)
    clsize(i) = sum(clust == i);
end
maxcsize = max(clsize);
return



function erase_string(str1)
fprintf(1,repmat('\b',1,length(str1))); % erase string
%fprintf(1,'\n');
return
