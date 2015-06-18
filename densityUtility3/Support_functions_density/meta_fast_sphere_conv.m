function indic = meta_fast_sphere_conv(xyzlist,xyzvox,r)
% indic = meta_fast_sphere_conv(xyzlist,xyzvox,r)
%
% convolution with a spherical filter
% output in indicator vector (vox x 1) for speed
% xyzlist is a list of voxel coords of in-mask voxels
% xyzvox is list of voxel coords for study (see meta_get_voxel_coords)
% r is radius in voxels
%
% see meta_read_mask to configure inputs

indic = zeros(size(xyzlist,1),1);

% get indicator of which voxels are in smoothed

for i = 1:size(xyzvox,1)
    
    [wh,indic]  = sphere_find_in_radius(xyzlist,xyzvox(i,:),5,indic);

end

return