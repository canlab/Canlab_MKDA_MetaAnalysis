function [mask,vol,voxsize] = get_random_mask(n,P)
%
% n is the number of random points
% P is the filename of the mask image file OR the actual mask array
%   P should be brain-only, or better, gray-matter only voxels
%   if P is empty, select using GUI
%
% vol is the search space - the whole brain mask where coordinates can be.
%
% structure of functions:
% 		voxel2mask.m			list of coordinates -> mask
%		mask2density.m			mask of coords -> density map
%
% for null hypothesis distribution generation:
%		generate_rdms.m	given: no. to generate, no. of points, mask size
%			get_random_mask.m	generate random list of coords and make into mask
%			mask2density.m			density map of null distribution
%			spm_write_vol.m		write output analyze image
%
%		density_npm_p.m			cluster p-value, and critical k 
%										given list of rdm images and threshold u			
%			count_clusters.m		density map -> no. of points, clusters 
%										at density threshold u
%										at a range of thresholds (input)
%										returns the number and avg. size of clusters
%			then finds p value based on distribution and 95% level

% ----------------------------------------------------
% * read the input brain mask image
% ----------------------------------------------------
if isempty(P)
    P = spm_get(1,'*.img','Select brain or gray matter mask');
end
if ischar(P)
    P = which(P);
    V = spm_vol(P);
    vol = spm_read_vols(V);
    voxsize = diag(V(1).mat)';
    voxsize = voxsize(1:3);
else
    vol = P;
    voxsize = NaN;
end



% ----------------------------------------------------
% * get list of all coordinates; select n at random
%   vol: row, col, array = x, y, z
% ----------------------------------------------------
[x,y,z] = ind2sub(size(vol),find(vol));     % find values > 0 in vol
allXYZ = [x y z];                           % XYZ is in 3 columns in this function
whichv = ceil(rand(1,n) * size(allXYZ,1));
XYZ = allXYZ(whichv,:);

ind = sub2ind(size(vol),XYZ(:,1),XYZ(:,2),XYZ(:,3));
mask = zeros(size(vol));
mask(ind) = 1;

% ----------------------------------------------------
% * deal with repeated coordinates by adding
% ----------------------------------------------------
ind = sort(ind);
repeats = ind(find(~diff(ind)));            % index values that are repeated
for i = 1:length(repeats)
    mask(repeats(i)) = mask(repeats(i)) + 1;
end

% this below may not allow mask values to be greater than 1 
% in the case of coordinate repeats.
% also probably too slow.  but no toolbox required!
%mask = voxel2mask(XYZ',size(vol));

return
