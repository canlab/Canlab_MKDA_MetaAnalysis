function [imgdata] = meta_read_image(name,mask)
% [imgdata] = meta_read_image(name,mask)
% 
% gets image data in a v x 1 voxel list (for fast processing)
% compatible with Meta_Activation_FWE subfunctions
%
% This is the inverse of meta_reconstruct_mask
% image (name) and mask (mask) must have same dimensions and origin!

if isempty(mask), mask = which('scalped_avg152T1_graymatter_smoothed.img');, end

V=spm_vol(mask); maskdata=spm_read_vols(V);

wh = find(maskdata);

V=spm_vol(name); imgdata=spm_read_vols(V);
imgdata = imgdata(wh);

return
