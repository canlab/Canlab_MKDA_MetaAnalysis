% [parcels, SVDinfo, parcel_stats] = Meta_Parcel(MC_Setup, varargin)
%
% documentation goes here.
%
%
% Example: use your own analysis mask:
% [parcels, SVDinfo, parcel_stats] = Meta_Parcel(MC_Setup, 'analysis_mask_name', 'acc_roi_mask.img');

function [parcels, SVDinfo, parcel_stats] = Meta_Parcel(MC_Setup, varargin)

global defaults
spm_defaults

%% user params: defaults

data_prctile = 90;          % save top xx % of voxels
n_components = 50;          % SVD components to use in parcellation
parcel_extent_thresh = 10;  % min # of contiguous voxels
analysis_mask_name = [];
dosave = 0;

mask = MC_Setup.volInfo.fname;
volInfo = MC_Setup.volInfo;

%% Optional inputs

for i = 1:length(varargin)
        if ischar(varargin{i})
            switch varargin{i}
              
                    % functional commands
                case 'data_prctile', data_prctile = varargin{i+1};
                case 'n_components', n_components = varargin{i+1};
                case 'parcel_extent_thresh', parcel_extent_thresh = varargin{i+1};
                case 'analysis_mask_name', analysis_mask_name = varargin{i+1};
                case {'dosave', 'save'}, dosave = 1;
                    
                otherwise, warning(['Unknown input string option:' varargin{i}]);
            end
        end
    end
%%

dat = MC_Setup.unweighted_study_data;

% get which voxels we care about
s = full(sum(dat'));

if ~isempty(analysis_mask_name)
    maskdat = scn_map_image(analysis_mask_name, mask);
    is_significant = find(maskdat(volInfo.wh_inmask));
else
    % Use top n%
    is_significant = find(s > prctile(s, data_prctile));
end

fprintf('Every voxel has at least %3.0f active contrasts\n', min(s(is_significant)));

fprintf('There are %3.0f voxels in the analysis area.\n', length(is_significant));


dat = dat(is_significant,:);
create_figure('# of Activating Contrasts'); plot(s); drawnow

%% singlar value D

dat = full(dat');
[U, S, V] = svd(dat, 'econ');
U = U(:,1:n_components); S = S(1:n_components, 1:n_components); V = V(:,1:n_components);
rec = U * S * V';


%% get parcels: max loading on SVs

nvox = size(dat, 2);
class = zeros(1, nvox);

for i = 1:nvox
        [taumtx] = correlation('tau', U, repmat(dat(:,i), 1, n_components));
        wh = find(taumtx == max(taumtx)); wh = wh(1);
        
        class(i) = wh;
        
        if mod(i, 100) == 0, fprintf('%3.0f ', i); end
        
end

%%
if dosave
    save SVD_data U S V dat class
end

%% Break into contiguous clusters

volInfo_sig = volInfo;
volInfo_sig.wh_inmask = volInfo_sig.wh_inmask(is_significant);
volInfo_sig.n_inmask = nvox;
volInfo_sig.xyzlist = volInfo_sig.xyzlist(is_significant, :);

%% Get parcels (cl structure) of contiguous voxels with same loading

parcels = {};
for i = 1:n_components
    myvec = class == i;
    parcels{i} = iimg_indx2clusters(myvec', volInfo_sig, .5, parcel_extent_thresh);

end

parcels = cat(2,parcels{:});

if dosave
    save SVD_data -append parcels volInfo_sig
end

%%

% get contrasts x 1 data indicator vector for each cluster 
% that specifies whether each contrast activated anywhere within cluster

[studybyparcel,studybyset] = Meta_cluster_tools('getdata',parcels, dat',volInfo_sig);

for i = 1:length(parcels), parcels(i).contrast_activation_dat = studybyparcel(:, i); end

[taumtx, t, p] = correlation('tau', studybyparcel);
p(p == 0) = 100*eps;

if dosave
    save SVD_data -append parcels studybyparcel studybyset taumtx t p
end

%% Threshold to get significance matrix - FDR

psq = p; psq(find(eye(size(p,1)))) = 0;
psq = squareform(psq);
pthr = FDR(p,.05);
if isempty(pthr), pthr = 0; end

sig = sign(t) .* (p < pthr);

sig(isnan(sig)) = 0;


%%
parcels(1).volInfo_sig = volInfo_sig;

SVDinfo = struct('U', U,  'S', S, 'V', V, 'class', class);
parcel_stats = struct('taumtx', taumtx, 't', t, 'p', p, 'sig', sig, 'fdr_threshold', pthr, 'studybyparcel', studybyparcel, 'studybyset', studybyset);

%% Make figures

try
    [parcel_stats, all_surf_handles] = Meta_Parcel_results(parcels, parcel_stats);
catch
    disp('Error displaying results');
end
    

end






