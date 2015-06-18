function indic = meta_stochastic_mask_blobs(volInfo,cl)
% indic = meta_stochastic_mask_blobs(volInfo,cl)
%
% Shuffles blob locations (centers) for a study
% Makes a mask first and then returns
% indicator (indic) format for aggregation across studies
%
% Inputs:
% volInfo.xyzlist: list of all in-mask coordinates (e.g.,from iimg_read_img)
% volInfo.dim: dimensions of mask 
% volInfo.nvox: number of voxels total in image
% volInfo.wh_inmask: indices of elements in mask (save only these)
%
% cl: clusters of contiguous voxels with XYZ field listing voxels
%
% Output:
% Vector of in-mask voxels (corresponding to xyzlist coordinates)
% With 1's where permuted blobs are.
% 
% tor wager


% ------------------------------------------------
% Set up input and bail out if no coordinates
% ------------------------------------------------

v = size(volInfo.xyzlist,1);
n = length(cl);

indic = zeros(volInfo.nvox,1);  % in full image space; reduce to in-mask later

if n == 0
    indic = indic(volInfo.wh_inmask);
    return
end

% ------------------------------------------------
% get a list of random blob centers
% ------------------------------------------------

wh = ceil(rand(n,1) .* v);
xyzvox = volInfo.xyzlist(wh,:);

% ------------------------------------------------
% add centers to each blob and add to mask
% ------------------------------------------------
xyzall = [];
for i = 1:n
    % add the randomly selected center to this blob
    xyz = cl(i).XYZ';
    
    if ~isempty(xyz)
        xyz(:,1) = xyz(:,1) + xyzvox(i,1);
        xyz(:,2) = xyz(:,2) + xyzvox(i,2);
        xyz(:,3) = xyz(:,3) + xyzvox(i,3);

        %emptymask(xyz(:,1),xyz(:,2),xyz(:,3)) = 1;  Soooooo slow
        xyzall = [xyzall; xyz];
    end

end

% ------------------------------------------------
% make sure all voxels are in range
% ------------------------------------------------

xyzall = min(xyzall, repmat(volInfo.dim(1:3),size(xyzall,1),1));
xyzall = max(xyzall, 1);

% ------------------------------------------------
% make index vector
% ------------------------------------------------

wh = sub2ind(volInfo.dim(1:3),xyzall(:,1),xyzall(:,2),xyzall(:,3));
indic(wh) = 1;

indic = indic(volInfo.wh_inmask);

return