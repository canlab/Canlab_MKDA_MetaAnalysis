function [maskdata,xyzlist,V] = meta_read_mask(mask)
% [maskdata,xyzlist,V] = meta_read_mask(mask)
% read mask and get XYZ voxel list of in-mask voxels
%
if isempty(mask), mask = which('scalped_avg152T1_graymatter_smoothed.img');, end

% make sure it exists
if ~exist(mask,'file'), disp(mask); error('Image does not exist!'); end

V=spm_vol(mask); maskdata=spm_read_vols(V);

% get voxel coordinates in mask
[x,y,z] = ind2sub(V.dim(1:3),find(maskdata));
xyzlist = [x y z];

return

% no masking -- slower
%
% dims = size(v);
% [x,y,z] = ind2sub(dims,(1:prod(dims))');
% xyzlist = [x y z];
% [wh,indic]  = sphere_find_in_radius(xyzlist,xyz,25);
% v2 = zeros(size(v));
% v2(wh) = 1;