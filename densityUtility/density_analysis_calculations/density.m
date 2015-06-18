function [dmt,clusters,dm,mask,CLU] = density(XYZ,u,varargin)
% function [dmt,clusters,dm,mask,CLU] = density(XYZ,u,[radius_mm],[output_filename],[zscores])
% 
% XYZ is 3-column vector of coordinates
%
% Tor Wager 2/18/02

P = ['brain_avg152T1.img'];     % 2 mm voxels
P = 'scalped_avg152T1_graymatter_smoothed.img';

if length(varargin) > 0, radius_mm = varargin{1}; 
else, radius_mm = 10;
end
disp(['Radius is ' num2str(radius_mm) ' mm.'])

if length(varargin) > 1, fname = varargin{2}; 
else, fname = input('Enter file name for the output image ','s');
end

if length(varargin) > 2,
    zscores = varargin{3};
    if length(zscores) ~= size(XYZ,1), error('Length of z-scores does not match number of points');,end
else
    zscores = 1;
end

if(length(fname) <= 4), error('Enter longer file name or type filename.img'),end
if ~strcmp(fname(end-4:end), '.img'), [d fname] = fileparts(fname);, fname = [fname '.img'], end

mask_file = P;
t1 = clock;

% -----------------------------------------------------
% * load standard brain
% -----------------------------------------------------
P = which(P);
V = spm_vol(P);
mask = zeros(V.dim(1:3));

voxsize = diag(V(1).mat)';
voxsize = voxsize(1:3);
radius = radius_mm ./ mean(voxsize);
sphere_vol = 4 * pi * radius_mm ^ 3 / 3;
XYZ = mm2vox(XYZ,V.mat);

% -----------------------------------------------------
% * make a mask out of XYZ, in space of P
% * deal with repeated coordinates by adding
% -----------------------------------------------------

ind = sub2ind(size(mask),XYZ(:,1),XYZ(:,2),XYZ(:,3));
mask(ind) = zscores;        % = 1 if no zscores entered, otherwise vector

ind = sort(ind);
repeats = ind(find(~diff(ind)));            % index values that are repeated
if length(zscores) > 1
    repeatz = zscores(find(~diff(ind)));            % z=scores of index values that are repeated
else
    repeatz = 1;
end
mask(repeats) = mask(repeats) + repeatz;

% -----------------------------------------------------
% * convert to density mask
% * write density mask before thresholding
% -----------------------------------------------------
dm = mask2density(mask,radius,[],sphere_vol);
%dm2 = mask2density(mask,radius);

%figure; subplot 211; imagesc(mask(:,:,33)); subplot 212; imagesc(dm(:,:,33))
%figure; subplot 211; imagesc(mask(:,:,33)); subplot 212; imagesc(dm(:,:,33).*sphere_vol)

V.fname = ['dens_' fname];
spm_write_vol(V,dm);

figure; hist(dm(dm>0),50)
if ~isempty(u),hold on; plot([u u],get(gca,'YLim'),'r','LineWidth',2),end

% -----------------------------------------------------
% * make and write filtered density mask
% -----------------------------------------------------
dmt = maskImg(dm,u,Inf);   
[d,fname,ext] = fileparts(fname);
V.fname = ['dens_' fname '_filtered' ext];
spm_write_vol(V,dmt);

% -----------------------------------------------------
% * get clusters from density mask for imaging
% -----------------------------------------------------
CLU = mask2struct(V.fname,u,0);
clusters = tor_extract_rois([],CLU,CLU);
eval(['save ' fname '_clusters CLU clusters'])

return



% -----------------------------------------------------
% * sub-functions
% -----------------------------------------------------

function XYZout = mm2vox(XYZ,M)
% converts a list of coordinates
% XYZ should be 3-column vector [x y z]
% calls tal2vox

XYZ = XYZ';
for i = 1:size(XYZ,2), 
	XYZout(i,:) = tal2vox(XYZ(:,i),M); 
end

XYZout = round(XYZout);

return

% -----------------------------------------------------
function vox=tal2vox(tal,M)
% converts from talairach coordinate to voxel coordinate
% based on variables from SPM.M (passed here for 
% faster operation)
% e.g., foo=tal2vox([-30 28 -30], VOL)
% from Russ Poldrack's spm ROI utility 

vox=[0 0 0];
vox(1)=(tal(1)-M(1,4))/M(1,1);
vox(2)=(tal(2)-M(2,4))/M(2,2);
vox(3)=(tal(3)-M(3,4))/M(3,3);

return





function dm = mask2density(mask,radius,varargin)
% function dm = mask2density(mask,radius,[opt] searchmask, [opt] sphere_vol)
%
% mask is the mask with ones where activation points are
% radius is in voxels
%
% optional arguments:
% 1 searchmask
%   searchmask is the whole brain search space [optional]
%   Mask and searchmask must be in same space, with same voxel sizes and dimensions!
%   If empty, search mask is ones(size(mask))
%
% 2 sphere_vol
%   sphere_vol is optional argument that translates count mask into density mask
%   by dividing the number of counts by this value.
%
% structure of functions:
% 		voxel2mask.m			list of coordinates -> mask
%		mask2density.m			mask of coords -> density map
%
% for null hypothesis distribution generation:
%		generate_rdms.m	given: no. to generate, no. of points, mask size
%	        get_random_mask.m	generate random list of coords
%			voxel2mask.m 
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

verbose = 0;
if verbose, t1 = clock;, end

if length(varargin) > 1
    sphere_vol = varargin{2};
else 
    sphere_vol = 0;
end

if length(varargin) > 0
    searchmask = varargin{1};
    if isempty(searchmask), searchmask = ones(size(mask));, end
else
    searchmask = ones(size(mask));
end

if any(size(mask) - size(searchmask))
    error('Mask and searchmask must be of the same dimensions and have the same voxel sizes!')
end

fmethod = 1;
switch fmethod
case 1
% --------------------------------------------------------------------------------------
% * define sphere to convolve with mask
% --------------------------------------------------------------------------------------
% MUCH slower.
sph = get_sphere(radius);
dm = convn(mask,sph,'same');

case 2

% --------------------------------------------------------------------------------------
% * define search space: coordinate list of whole mask
% --------------------------------------------------------------------------------------
if verbose, fprintf(1,'defining search space ... '), end
% [x,y,z] = ind2sub(size(searchmask),find(searchmask));     % find coords of values > 0 in mask
[x,y,z] = ind2sub(size(mask),find(searchmask > 0));     % find coords of values > 0 in mask
vindex = [x y z]';

% old way: slower
% vindex  = build_voxel_index(mask);

Q       = ones(1,size(vindex,2));
dm      = zeros(size(mask));
slices  = vindex(3,:);                              % used to pass only certain slices to subfunctions


% --------------------------------------------------------------------------------------
% * locate coordinates in mask to make spheres around
% --------------------------------------------------------------------------------------
[x,y,z] = ind2sub(size(mask),find(mask));     % find coordinates of values > 0 in mask
vr = [x y z]';

if verbose, fprintf(1,'\nrestricting to spheres around %3.0f points ... ',size(vr,2)), end

% the idea here is to pass into add_density_sphere only coordinates
% that are in slices that could possibly be within the sphere, to speed things up.
% this next line tells us when we need to update to new slices
% when the slice changes, we need a new set of slices to pass to add_density_sphere
% vr should already be sorted ascending wrt slices, so no need to sort again.
upd = [1 diff(vr(3,:)) > 0];

npts = size(vr,2);
str = fprintf(1,'Computing %3.0f points: %03d',npts,0);

for i = 1:npts                                               % for each reported point...
    
    fprintf(1,'\b\b\b%03d',i);
    
    if upd(i)
        vin = vindex(:,slices <= vr(3,i) + radius & slices > vr(3,i) - radius);    % these are coordinates that may be in-sphere
    end
    
    count   = mask(vr(1,i),vr(2,i),vr(3,i)); 
										% get count at this point
    dm      = add_density_sphere(vr(:,i),vin,Q,count,dm,radius);        % add density to sphere around this point
end

erase_string(str);

% --------------------------------------------------------------------------------------
% * turn into density -> count per unit area of sphere 
% --------------------------------------------------------------------------------------
if sphere_vol
    dm = dm ./ sphere_vol;
end

if verbose, 
    t2 = clock;
    fprintf(1,'\ndone in %3.0f s!',etime(t2,t1))
end

end

return





% Sub-functions
% --------------------------------------------------------------------------------------
function vindex = build_voxel_index(mask)
% NOT USED
rows = size(mask,1);
cols = size(mask,2);
a = []; b = []; c = []; vindex = [];
for i = 1:cols
        a = [a 1:rows];
        b = [b i * ones(1,rows)];
end
c = [a;b];
for i = 1:size(mask,3)
    c(3,:) = i;
    vindex = [vindex c];
end
return


function [dens,j] = get_density(coord,vindex,Q,mask,radius,sphere_vol)
    % NOT USED
    % j is index numbers of in-sphere coordinates
    % this line is slow for big vindexes (1 mil).  4.17 s
    j       = find(sum((vindex - coord*Q).^2) <= radius);

    %	total   = length(j) * voxel_vol;
    count   = sum(mask(j) > 0);
    dens    = count / sphere_vol;
return


function [dm,j] = add_density_sphere(coord,vindex,Q,count,dm,radius)
    % j is index numbers of in-sphere coordinates
    % coord is coordinate of sphere center
    % vindex is XYZ index of all slices in the vicinity
    % Q is ones of (3,length vindex)
    % count is the number of counts at the sphere center
    
    j       = find(sum((vindex - coord*Q(:,1:size(vindex,2))).^2) <= radius^2);

    % turn index of restricted vindex space back to whole-mask coordinates
    XYZ = vindex(:,j);
    ind = sub2ind(size(dm),XYZ(1,:),XYZ(2,:),XYZ(3,:));
    
    % add the mask value at the point to voxels in sphere; mask may be > 1
    % do not divide by sphere_vol - instead, do this at the end.
    dm(ind) = dm(ind) + count;
    
    
return


function sph_mask = get_sphere(radius)
    % returns spherical mask in 3-d array, with unit density in elements
    % NOT USED
    sph_mask    = zeros(2*radius+1,2*radius+1,2*radius+1);
    vindex      = build_voxel_index(sph_mask);
    Q           = ones(1,size(vindex,2));
    sph_mask(radius,radius,radius) = 1;
    [dens,j]    = get_density([radius;radius;radius],vindex,Q,sph_mask,radius,1);
    sph_mask(j) = dens;
return


function erase_string(str)
    fprintf(1,repmat('\b',1,length(str))); % erase string
    %fprintf(1,'\n');
return
